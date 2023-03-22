from copy import deepcopy
from pathlib import Path
import numpy as np
from astropy.table import vstack
import astropy.units as u
from astropy.io import fits
import logging

from pyirf.binning import (
    create_bins_per_decade,
    add_overflow_bins,
)

from pyirf.io import create_aeff2d_hdu
from pyirf.utils import calculate_theta

from pyirf.irf import effective_area_per_energy

from pyirf.benchmarks import energy_bias_resolution, angular_resolution
from pyirf.benchmarks.energy_bias_resolution import energy_resolution_absolute_68

from lstchain.reco.utils import reco_source_position_sky
from lstchain.io import read_mc_dl2_to_QTable

import operator
from lstchain.io.event_selection import DL3Cuts, calculate_percentile_cut, evaluate_binned_cut

from ctapipe.core.traits import Float, Unicode


MIN_NEVENTS_RES = 100  # min number of events to consider the resolution in the plot
MIN_ENERGY = 20.0 * u.GeV
MAX_ENERGY = 20.05 * u.TeV
N_BIN_PER_DECADE = 5
TRUE_ENERGY_BINS = add_overflow_bins(create_bins_per_decade(MIN_ENERGY, MAX_ENERGY, N_BIN_PER_DECADE)).to(u.TeV)

logging.basicConfig(filename='dl2_to_irfs.log', filemode='a', format='%(name)s - %(levelname)s - %(message)s')



def patch_events_altaz(events, effective_focal_length=29.04 * u.m):
    """
    Recompute alt/az with the effective focal length
    """
    for dist in ['x', 'y', 'reco_disp_dx', 'reco_disp_dy']:
        if not events[dist].unit:
            events[dist].unit = u.m
    for ang in ['alt_tel', 'az_tel']:
        if not events[ang].unit:
            events[ang].unit = u.rad

    raltaz = reco_source_position_sky(events['x'],
                                      events['y'],
                                      events['reco_disp_dx'],
                                      events['reco_disp_dy'],
                                      effective_focal_length,
                                      events['alt_tel'],
                                      events['az_tel'])

    events['reco_alt'] = raltaz.alt.to(u.rad)
    events['reco_az'] = raltaz.az.to(u.rad)
    events['gammaness'] = events['gh_score']
    events["theta"] = calculate_theta(events, 
                                      assumed_source_az=events['true_az'],
                                      assumed_source_alt=events['true_alt'],
                                     )

    return None


def load_event_dict(filename, events_dict, force_reload=False):
    if not isinstance(events_dict, dict):
        raise TypeError("events_dict must be a dict")
    node = Path(filename).parent.name
    if node not in events_dict or force_reload:
        print(f"loading node {node}")
        events_dict[node] = {'path': Path(filename)}
        events_dict[node]['params'], events_dict[node]['simu_info'] = read_mc_dl2_to_QTable(filename)
        # Patch for the coma aberration bias
        patch_events_altaz(events_dict[node]['params'])
    else:
        print("loaded already")


def groupby_alt(pointings):
    pp = pointings.to_pandas()
    return pp.groupby('alt').groups


def groupby_az_sym_magnorth(pointings):
    """
    group pointings by azimuths symmetrical to magnetic north
    """
    pp = pointings.to_pandas()
    pp['cosaz'] = np.cos(np.deg2rad(pp['az'] + 180 - 175.158)).round(3)
    return pp.groupby('cosaz').groups


def groupby_alt_az_sym_magnorth(pointings):
    """
    group pointings by azimuths symmetrical to magnetic north
    """
    pp = pointings.to_pandas()
    pp['cosaz'] = np.cos(np.deg2rad(pp['az'] + 180 - 175.158)).round(3)
    return pp.groupby('cosaz').groups


def merge_per_alt(events_dict, pointings):
    stacked_alt = {}
    grp = groupby_alt(pointings)

    for alt, index in grp.items():
        print(alt)
        stacked_alt[alt] = {}
        stacked_alt[alt]['params'] = vstack([events_dict[p['dirname']]['params'] for p in pointings[list(index)]])
        stacked_alt[alt]['simu_info'] = deepcopy(events_dict[pointings[list(index)[0]]['dirname']]['simu_info'])
        stacked_alt[alt]['simu_info'].n_showers = sum(
            [events_dict[p['dirname']]['simu_info'].n_showers for p in pointings[list(index)]])

    return stacked_alt


def calculate_angular_resolution(gammas, true_energy_bins=TRUE_ENERGY_BINS):
    ang_res = angular_resolution(gammas, true_energy_bins)
    n, bins = np.histogram(gammas['true_energy'], bins=true_energy_bins)
    ang_res['angular_resolution'][n < MIN_NEVENTS_RES] = np.nan
    return ang_res


def calculate_energy_bias_resolution(gammas, true_energy_bins=TRUE_ENERGY_BINS):
    bias_resolution = energy_bias_resolution(gammas,
                                             true_energy_bins,
                                             resolution_function=energy_resolution_absolute_68)

    n, bins = np.histogram(gammas['true_energy'], bins=true_energy_bins)
    bias_resolution['resolution'][n < MIN_NEVENTS_RES] = np.nan
    return bias_resolution


def filter_intensity(events, intensity=50):
    return events[events['intensity'] > intensity]


def filter_wrong_tail(events):
    return events[events['reco_disp_sign'] == events['disp_sign']]


class DL3CutsIRFs(DL3Cuts):
    """
    redefine some methods to apply cuts on true energy
    """

    alpha_containment = Float(
        help="Percentage containment region for alpha cuts",
        default=1.0,
    ).tag(config=True)

    min_alpha_cut = Float(
        help="Minimum alpha cut (deg) in an energy bin",
        default_value=0,
    ).tag(config=True)

    max_alpha_cut = Float(
        help="Maximum alpha cut (deg) in an energy bin",
        default_value=90,
    ).tag(config=True)

    fill_alpha_cut = Float(
        help="Fill value of alpha cut (deg) in an energy bin with fewer " +
             "than minimum number of events present",
        default_value=90,
    ).tag(config=True)
    
    energy_type = Unicode(
        help="reco or true",
        default_value="true",
    ).tag(config=True)
    

    def energy_dependent_gh_cuts(self, data, energy_bins, smoothing=None):
        """
        Evaluating energy-dependent gammaness cuts, in a given
        data, with provided reco energy bins, and other parameters to
        pass to the pyirf.cuts.calculate_percentile_cut function
        """

        gh_cuts = calculate_percentile_cut(
            data["gh_score"],
            data[f"{self.energy_type}_energy"],
            bins=energy_bins,
            min_value=self.min_gh_cut,
            max_value=self.max_gh_cut,
            fill_value=data["gh_score"].max(),
            percentile=100 * (1 - self.gh_efficiency),
            smoothing=smoothing,
            min_events=self.min_event_p_en_bin,
        )
        return gh_cuts

    def apply_energy_dependent_gh_cuts(self, data, gh_cuts):
        """
        Applying a given energy-dependent gh cuts to a data file, along the reco
        energy bins provided.
        """

        data["selected_gh"] = evaluate_binned_cut(
            data["gh_score"],
            data[f"{self.energy_type}_energy"],
            gh_cuts,
            operator.ge,
        )
        return data[data["selected_gh"]]

    def energy_dependent_theta_cuts(
            self, data, energy_bins, smoothing=None,
    ):
        """
        Evaluating an optimized energy-dependent theta cuts, in a given
        data, with provided reco energy bins, and other parameters to
        pass to the pyirf.cuts.calculate_percentile_cut function.
        Note: Using too fine binning will result in too un-smooth cuts.
        """

        theta_cuts = calculate_percentile_cut(
            data["theta"],
            data[f"{self.energy_type}_energy"],
            bins=energy_bins,
            min_value=self.min_theta_cut * u.deg,
            max_value=self.max_theta_cut * u.deg,
            fill_value=self.fill_theta_cut * u.deg,
            percentile=100 * self.theta_containment,
            smoothing=smoothing,
            min_events=self.min_event_p_en_bin,
        )
        return theta_cuts

    def apply_energy_dependent_theta_cuts(self, data, theta_cuts):
        """
        Applying a given energy-dependent theta cuts to a data file, along the
        reco energy bins provided.
        """

        data["selected_theta"] = evaluate_binned_cut(
            data["theta"],
            data[f"{self.energy_type}_energy"],
            theta_cuts,
            operator.le,
        )
        return data[data["selected_theta"]]

    def energy_dependent_alpha_cuts(
            self, data, energy_bins, smoothing=None
    ):
        """
        Evaluating an optimized energy-dependent alpha cuts, in a given
        data, with provided reco energy bins, and other parameters to
        pass to the pyirf.cuts.calculate_percentile_cut function.
        Note: Using too fine binning will result in too un-smooth cuts.
        """

        alpha_cuts = calculate_percentile_cut(
            data["alpha"],
            data[f"{self.energy_type}_energy"],
            bins=energy_bins,
            min_value=self.min_alpha_cut * u.deg,
            max_value=self.max_alpha_cut * u.deg,
            fill_value=self.fill_alpha_cut * u.deg,
            percentile=100 * self.alpha_containment,
            smoothing=smoothing,
            min_events=self.min_event_p_en_bin,
        )
        return alpha_cuts

    def apply_energy_dependent_alpha_cuts(self, data, alpha_cuts):
        """
        Applying a given energy-dependent alpha cuts to a data file, along the
        reco energy bins provided.
        """

        data["selected_alpha"] = evaluate_binned_cut(
            data["alpha"],
            data[f"{self.energy_type}_energy"],
            alpha_cuts,
            operator.le,
        )
        return data[data["selected_alpha"]]


def filter_events(events,
                  gh_efficiency=1.0,
                  min_theta_cut=0.0 * u.deg,
                  max_theta_cut=np.inf * u.deg,
                  global_theta_cut=np.inf * u.deg,
                  global_alpha_cut=np.inf * u.deg,
                  min_gh_cut=0.,
                  max_gh_cut=1.0,
                  global_gh_cut=1.0,
                  min_event_p_en_bin=100,
                  theta_containment=1.0,
                  alpha_containment=1.0,
                  srcdep=False,
                  filter_disp_tail=False,
                  energy_type='true',
                  ):
    events = filter_intensity(events)
    
    if srcdep and filter_disp_tail:
        raise ValueError("You don't want to filter for wrong tail assignment in case of source-dep analysis")


    dl3_cuts = DL3CutsIRFs(
        min_event_p_en_bin=min_event_p_en_bin,
        min_gh_cut=min_gh_cut,
        max_gh_cut=max_gh_cut,
        global_gh_cut=global_gh_cut,
        gh_efficiency=gh_efficiency,
        min_theta_cut=min_theta_cut.to_value(u.deg),
        max_theta_cut=max_theta_cut.to_value(u.deg),
        fill_theta_cut=max_theta_cut.to_value(u.deg),
        global_theta_cut=global_theta_cut.to_value(u.deg),
        theta_containment=theta_containment,
        min_alpha_cut=0,
        max_alpha_cut=90,
        global_alpha_cut=global_alpha_cut.to_value(u.deg),
        fill_alpha_cut=global_alpha_cut.to_value(u.deg),
        alpha_containment=alpha_containment,
        energy_type=energy_type,
    )

    gh_cuts = dl3_cuts.energy_dependent_gh_cuts(events, energy_bins=TRUE_ENERGY_BINS)
    events = dl3_cuts.apply_energy_dependent_gh_cuts(events, gh_cuts)

    if srcdep:
        alpha_cuts = dl3_cuts.energy_dependent_alpha_cuts(events, energy_bins=TRUE_ENERGY_BINS)
        events = dl3_cuts.apply_energy_dependent_alpha_cuts(events, alpha_cuts)
        return events, gh_cuts, alpha_cuts

    else:
        filtered_wrong_tail_events = filter_wrong_tail(events)
        # theta_cuts = dl3_cuts.energy_dependent_theta_cuts(events, energy_bins=TRUE_ENERGY_BINS)
        theta_cuts = dl3_cuts.energy_dependent_theta_cuts(filtered_wrong_tail_events, energy_bins=TRUE_ENERGY_BINS)
        
        if filter_disp_tail:
            # we filter for wrong disp only for AngRes, where we don't apply theta cuts
            logging.warning("filtering wrong disp sign events")
            events = filtered_wrong_tail_events
        else:
            # we apply theta cuts only for Aeff and Eres, where we don't filter for wrong disp (but the theta cuts are still computed without taking these wrongly assigned events)
            logging.warning("applying theta cuts")
            events = dl3_cuts.apply_energy_dependent_theta_cuts(events, theta_cuts)
    return events, gh_cuts, theta_cuts


def produce_irfs_pyirf(gammas, simu_info, gh_efficiency, outfile,
                       theta_containment=1.0,
                       alpha_containment=1.0,
                       srcdep=False,
                       energy_type='reco'
                       ):


    filter_disp_tail = True if not srcdep else False
    gammas_ang_res, gh_cuts, theta_cuts_ang_res = filter_events(gammas,
                                                                gh_efficiency=gh_efficiency,
                                                                srcdep=srcdep,
                                                                filter_disp_tail=filter_disp_tail,
                                                                energy_type=energy_type,
                                                                )

    gammas_energy_effarea, gh_cuts_, theta_cuts = filter_events(gammas,
                                                                gh_efficiency=gh_efficiency,
                                                                theta_containment=theta_containment,
                                                                alpha_containment=alpha_containment,
                                                                srcdep=srcdep,
                                                                filter_disp_tail=False,
                                                                energy_type=energy_type,
                                                                )

    logging.warning(f"{len(gammas)} gammas --> {len(gammas_ang_res)} ({len(gammas_ang_res)/len(gammas):.2f}) --> {len(gammas_energy_effarea)} ({len(gammas_energy_effarea)/len(gammas):.2f})")

    gh_cuts.rename_columns(['low', 'high', 'center', 'cut'],
                           ['true_energy_low', 'true_energy_high', 'true_energy_center', 'gh_cuts'])

    gh_table = gh_cuts

    hdus = [
        fits.PrimaryHDU(),
        fits.BinTableHDU(gh_table, name="GH_CUTS"),
        fits.BinTableHDU(theta_cuts, name="THETA_CUTS"),
        fits.BinTableHDU(theta_cuts_ang_res, name="THETA_CUTS_ANG_RES"),
    ]

    masks = {
        "_NO_CUTS": gammas,
        "_ANG_RES": gammas_ang_res,
        "": gammas_energy_effarea,
    }

    for label, selected_gammas in masks.items():
        effective_area = effective_area_per_energy(
            selected_gammas,
            simu_info,
            true_energy_bins=TRUE_ENERGY_BINS,
        )
        hdus.append(
            create_aeff2d_hdu(
                effective_area[..., np.newaxis],  # add one dimension for FOV offset
                TRUE_ENERGY_BINS,
                [0, 10] * u.deg,
                extname="EFFECTIVE_AREA" + label,
            )
        )
    ang_res = calculate_angular_resolution(gammas_ang_res, TRUE_ENERGY_BINS)
    bias_resolution = calculate_energy_bias_resolution(gammas_energy_effarea, TRUE_ENERGY_BINS)

    hdus.append(fits.BinTableHDU(ang_res, name="ANGULAR_RESOLUTION"))
    hdus.append(fits.BinTableHDU(bias_resolution, name="ENERGY_BIAS_RESOLUTION"))

    Path(outfile).parent.mkdir(exist_ok=True)
    fits.HDUList(hdus).writeto(outfile, overwrite=True)

