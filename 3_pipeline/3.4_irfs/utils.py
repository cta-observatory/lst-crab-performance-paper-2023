import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import ctaplot
from ctaplot.ana import logbin_mean
from astropy.visualization import quantity_support

from pyirf.binning import (
    create_bins_per_decade,
    add_overflow_bins,
)
from pyirf.benchmarks import energy_bias_resolution, angular_resolution
from pyirf.benchmarks.energy_bias_resolution import energy_resolution_absolute_68
from pyirf.utils import calculate_theta
from lstchain.reco.utils import reco_source_position_sky
from pyirf import irf
from lstchain.io import read_mc_dl2_to_QTable

from pathlib import Path
from astropy.table import QTable


MIN_THETA_CUT = 0.1 * u.deg
MAX_THETA_CUT = 0.5 * u.deg

MIN_ENERGY = 20.0 * u.GeV
MAX_ENERGY = 20.05 * u.TeV

N_BIN_PER_DECADE = 5

MIN_NEVENTS_RES = 100  #min number of events to consider the resolution in the plot

default_res_plot_opt = {
    'elinewidth': 1,
    'ls': '-',
    'markersize': 1,
    'fmt': 'o',
}


# binnings for the irfs
true_energy_bins = add_overflow_bins(create_bins_per_decade(MIN_ENERGY, MAX_ENERGY, N_BIN_PER_DECADE)).to(u.TeV)
reco_energy_bins = add_overflow_bins(create_bins_per_decade(MIN_ENERGY, MAX_ENERGY, N_BIN_PER_DECADE)).to(u.TeV)
true_energy_bins_center = np.sqrt(true_energy_bins[:-1]*true_energy_bins[1:])
# true_energy_bins_center = logbin_mean(true_energy_bins)
true_energy_bins_xerr = (true_energy_bins_center - true_energy_bins[:-1], true_energy_bins[1:] - true_energy_bins_center)

# true_energy_bins = ctaplot.irf_cta().energy_bins
# reco_energy_bins = true_energy_bins


# fov_offset_bins = [0, 0.6] * u.deg
# source_offset_bins = np.arange(0, 1 + 1e-4, 1e-3) * u.deg
# energy_migration_bins = np.geomspace(0.2, 5, 200)




def patch_events_altaz(events, effective_focal_length=29.04*u.m):
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
    # events['offset'] = ctaplot.ana.angular_separation_altaz(events['reco_alt'],
    #                                                         events['reco_az'],
    #                                                         events['true_alt'],
    #                                                         events['true_az'])
    events["theta"] = calculate_theta(events, assumed_source_az=events['true_az'], assumed_source_alt=events['true_alt'])

    return None



def efficiency_to_gammaness(gammaness, efficiency):
    sorted_gam = np.sort(gammaness)
    index = int((1-efficiency) * len(sorted_gam))
    return sorted_gam[index]


def filter_efficiency(events, efficiency):
    gam_cut = efficiency_to_gammaness(events['gammaness'], efficiency)
    return events[events['gammaness']>gam_cut]


def plot_effective_area(gammas, simu_info, ax=None, **kwargs):
    ax = plt.gca() if ax is None else ax
    effective_area = irf.effective_area_per_energy(
        gammas,
        simu_info,
        true_energy_bins=true_energy_bins,
    )
    with quantity_support():
        ax.stairs(effective_area,  edges=true_energy_bins, **kwargs)
    ax.set_ylabel(r'Effective Area [$m^2$]')
    ax.set_xlabel(r'$E_T$ [TeV]')
    # ax.grid(True, which='both')
    # ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Effective area')
    return ax


def plot_energy_resolution(gammas, ax=None, **kwargs):
    ax = plt.gca() if ax is None else ax
    bias_resolution = energy_bias_resolution(gammas,
                                             true_energy_bins,
                                             resolution_function=energy_resolution_absolute_68)
    
    n, bins = np.histogram(gammas['true_energy'], bins=true_energy_bins)
    bias_resolution['resolution'][n < MIN_NEVENTS_RES] = np.nan
    
    for k, v in default_res_plot_opt.items():
         kwargs.setdefault(k, v)
            
    with quantity_support():
        ax.errorbar(true_energy_bins_center, bias_resolution['resolution'], xerr=true_energy_bins_xerr, **kwargs)
    # ax.legend()
    ax.set_xscale('log')
    ax.grid(True, which='both')
    ax.set_ylim(0, 0.6)
    ax.set_title('Energy resolution')
    ax.set_ylabel(r'$\left(\Delta E/E\right)_{68}$')
    ax.set_xlabel(r'$E_T$ [TeV]')
    return ax


def plot_energy_bias(gammas, ax=None, **kwargs):
    ax = plt.gca() if ax is None else ax
    bias_resolution = energy_bias_resolution(gammas,
                                             true_energy_bins,
                                             resolution_function=energy_resolution_absolute_68)
    n, bins = np.histogram(gammas['true_energy'], bins=true_energy_bins)
    bias_resolution['bias'][n < MIN_NEVENTS_RES] = np.nan
    
   
    for k, v in default_res_plot_opt.items():
         kwargs.setdefault(k, v)
    with quantity_support():
        ax.errorbar(true_energy_bins_center, bias_resolution['bias'], xerr=true_energy_bins_xerr, **kwargs)
    # ax.legend()
    ax.set_xscale('log')
    ax.grid(True, which='both')
    ax.set_title('Energy bias')
    ax.set_ylabel('bias')
    ax.set_xlabel(r'$E_T$ [TeV]')
    return ax



def calculate_angular_resolution(gammas):
    ang_res = angular_resolution(gammas, true_energy_bins)
    n, bins = np.histogram(gammas['true_energy'], bins=true_energy_bins)
    ang_res['angular_resolution'][n < MIN_NEVENTS_RES] = np.nan
    return ang_res

def save_angular_resolution(ang_res):
    ang_res


def plot_angular_resolution(gammas, ax=None, **kwargs):
    ax = plt.gca() if ax is None else ax
    ang_res = angular_resolution(gammas, true_energy_bins)
    n, bins = np.histogram(gammas['true_energy'], bins=true_energy_bins)
    ang_res['angular_resolution'][n < MIN_NEVENTS_RES] = np.nan
    
    for k, v in default_res_plot_opt.items():
         kwargs.setdefault(k, v)
            
    with quantity_support():
        ax.errorbar(true_energy_bins_center, ang_res['angular_resolution'], xerr=true_energy_bins_xerr, **kwargs)
    # ax.legend()
    ax.set_xscale('log')
    ax.grid(True, which='both')
    ax.set_title('Angular resolution')
    ax.set_ylabel('Angular resolution [degrees]')
    ax.set_xlabel(r'$E_T$ [TeV]')
    return ax


def find_closests_test_pointings(dec='dec_2276'):
    from lstmcpipe.config import paths_config
    pc = paths_config.PathConfigAllSkyFull('test', [dec])
    
    closests = []
    for point in pc.train_configs[dec].pointings:
        ii = np.argmin(ctaplot.ana.angular_separation_altaz(point['alt'], point['az'], pc.test_configs[dec].pointings['alt'], pc.test_configs[dec].pointings['az']))
        closests.append([pc.test_configs[dec].pointings[ii]['alt'].to_value(u.rad), pc.test_configs[dec].pointings[ii]['az'].to_value(u.rad),
                         pc.test_configs[dec].pointings[ii]['dirname']
                        ])
    closests, unique_index = np.unique(closests, axis=0, return_index=True)
    return closests


def filter_theta_cut(events, theta_cut):
    return events[events['theta']<theta_cut]

def filter_efficiency_true_energy_dep(events, efficiency, true_energy_bins):
    indices = np.digitize(events['true_energy'], true_energy_bins)
    mask = np.zeros(len(events), dtype=bool)
    gammaness_cuts = []
    for ii in range(1, len(true_energy_bins)):
        gam_cut = efficiency_to_gammaness(events[indices == ii]['gammaness'], efficiency)
        mask[(indices == ii) & (events['gammaness'] > gam_cut)] = True
        gammaness_cuts.append(gam_cut)
    return events[mask], gammaness_cuts
        
def filter_intensity(events, intensity=50):
    return events[events['intensity']>intensity]


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
    pp['cosaz'] = np.cos(np.deg2rad(pp['az']+180-175.158)).round(3)
    return pp.groupby('cosaz').groups


def groupby_alt_az_sym_magnorth(pointings):
    """
    group pointings by azimuths symmetrical to magnetic north
    """
    pp = pointings.to_pandas()
    pp['cosaz'] = np.cos(np.deg2rad(pp['az']+180-175.158)).round(3)
    return pp.groupby('cosaz').groups


def plot_gh_cut_per_energy(filename, ax=None, **kwargs):
    """ """
    ax = plt.gca() if ax is None else ax

    gh_cut = QTable.read(filename, hdu="GH_CUTS")[1:-1]

    kwargs.setdefault("ls", "")
    ax.errorbar(
        0.5 * (gh_cut["true_energy_low"] + gh_cut["true_energy_high"]).to_value(u.TeV),
        gh_cut["gh_cuts"],
        xerr=0.5 * (gh_cut["true_energy_high"] - gh_cut["true_energy_low"]).to_value(u.TeV),
        **kwargs,
    )

    ax.legend()
    ax.set_ylabel("G/H-cut")
    ax.set_xlabel(r"$E_\mathrm{reco} / \mathrm{TeV}$")
    ax.set_xscale("log")

    return 


def read_angular_resolution(filename):
    ang_res = QTable.read(filename, hdu="ANGULAR_RESOLUTION")[1:-1]
    ang_res.round(3)
    return ang_res


def read_energy_resolution_bias(filename):
    ene_res = QTable.read(filename, hdu="ENERGY_BIAS_RESOLUTION")[1:-1]
    ene_res.round(3)
    return ene_res


def read_effective_area(filename):
    eff_area = QTable.read(filename, hdu="EFFECTIVE_AREA")[0]
    return eff_area


def plot_theta_cut_per_energy(filename, hdu="THETA_CUTS", ax=None, **kwargs):
    """ """
    ax = plt.gca() if ax is None else ax

    th_cut = QTable.read(filename, hdu=hdu)[1:-1]

    kwargs.setdefault("ls", "")
    ax.errorbar(
        0.5 * (th_cut["low"] + th_cut["high"]).to_value(u.TeV),
        th_cut["cut"],
        xerr=0.5 * (th_cut["high"] - th_cut["low"]).to_value(u.TeV),
        **kwargs,
    )

    ax.legend()
    ax.set_ylabel(r"$\theta$ cut")
    ax.set_xlabel(r"$E_\mathrm{reco} / \mathrm{TeV}$")
    ax.set_xscale("log")

    return