#!/usr/bin/env python
from astropy.table import QTable
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u


def plot_gh_cut_per_energy(filename, ax=None, **kwargs):
    """ """
    ax = plt.gca() if ax is None else ax

    gh_cut = QTable.read(filename, hdu="GH_CUTS")[1:-1]

    kwargs.setdefault("ls", "")
    ax.errorbar(
        0.5 * (gh_cut["low"] + gh_cut["high"]).to_value(u.TeV),
        gh_cut["cut"],
        xerr=0.5 * (gh_cut["high"] - gh_cut["low"]).to_value(u.TeV),
        **kwargs,
    )

    ax.set_ylabel("G/H-cut")
    ax.set_xlabel(r"$E_\mathrm{reco} / \mathrm{TeV}$")
    ax.set_xscale("log")
    ax.grid(True)

    return ax



def plot_effective_area_from_file(file, all_cuts=False, ax=None, **kwargs):
    """ """

    ax = plt.gca() if ax is None else ax

    if all_cuts:
        names = ["", "_NO_CUTS", "_ONLY_GH", "_ONLY_THETA"]
    else:
        names = tuple([""])

    label_basename = kwargs["label"] if "label" in kwargs else ""

    kwargs.setdefault("ls", "")

    for name in names:

        area = QTable.read(file, hdu="EFFECTIVE_AREA" + name)[0]

        kwargs["label"] = label_basename + name.replace("_", " ")

        ax.errorbar(
            0.5 * (area["ENERG_LO"] + area["ENERG_HI"]).to_value(u.TeV)[1:-1],
            area["EFFAREA"].to_value(u.m**2).T[1:-1, 0],
            xerr=0.5 * (area["ENERG_LO"] - area["ENERG_HI"]).to_value(u.TeV)[1:-1],
            **kwargs,
        )

        # Style settings
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$E_\mathrm{True}$ / TeV")
        ax.set_ylabel("Effective collection area / m²")
        ax.grid(which="both")
        ax.grid(True, which="both")

    return ax


def plot_psf_from_file(filename):

    psf_table = QTable.read(filename, hdu="PSF")[0]
    # select the only fov offset bin
    psf = psf_table["RPSF"].T[:, 0, :].to_value(1 / u.sr)

    offset_bins = np.append(psf_table["RAD_LO"], psf_table["RAD_HI"][-1])
    phi_bins = np.linspace(0, 2 * np.pi, 100)

    # Let's make a nice 2d representation of the radially symmetric PSF
    r, phi = np.meshgrid(offset_bins.to_value(u.deg), phi_bins)

    # look at a single energy bin
    # repeat values for each phi bin
    center = 0.5 * (psf_table["ENERG_LO"] + psf_table["ENERG_HI"])

    fig = plt.figure(figsize=(15, 5))
    axs = [fig.add_subplot(1, 3, i, projection="polar") for i in range(1, 4)]

    for bin_id, ax in zip([10, 20, 30], axs):
        image = np.tile(psf[bin_id], (len(phi_bins) - 1, 1))

        ax.set_title(f"PSF @ {center[bin_id]:.2f} TeV")
        ax.pcolormesh(phi, r, image)
        ax.set_ylim(0, 0.25)
        ax.set_aspect(1)

    fig.tight_layout()

    return axs


def plot_theta_cut_from_file(filename, ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax

    rad_max = QTable.read(filename, hdu="RAD_MAX")[0]

    kwargs.setdefault("ls", "")
    ax.errorbar(
        0.5 * (rad_max["ENERG_LO"] + rad_max["ENERG_HI"])[1:-1].to_value(u.TeV),
        rad_max["RAD_MAX"].T[1:-1, 0].to_value(u.deg),
        xerr=0.5 * (rad_max["ENERG_HI"] - rad_max["ENERG_LO"])[1:-1].to_value(u.TeV),
        **kwargs,
    )

    ax.set_ylabel("θ-cut / deg²")
    ax.set_xlabel(r"$E_\mathrm{reco} / \mathrm{TeV}$")
    ax.set_xscale("log")

    return ax


def plot_angular_resolution_from_file(filename, ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax

    ang_res = QTable.read(filename, hdu="ANGULAR_RESOLUTION")[1:-1]

    kwargs.setdefault("ls", "")
    ax.errorbar(
        0.5 * (ang_res["true_energy_low"] + ang_res["true_energy_high"]).to_value(u.TeV),
        ang_res["angular_resolution"].to_value(u.deg),
        xerr=0.5 * (ang_res["true_energy_high"] - ang_res["true_energy_low"]).to_value(u.TeV),
        **kwargs,
    )

    # Style settings
    #     ax.set_xlim(1.e-2, 2.e2)
    #     ax.set_ylim(2.e-2, 1)
    ax.set_xscale("log")
    ax.set_xlabel(r"$E_\mathrm{True}$ / TeV")
    ax.set_ylabel("Angular Resolution / deg")
    ax.grid(True, which="both")

    return ax


def plot_energy_dispersion_from_file(filename):

    edisp = QTable.read(filename, hdu="ENERGY_DISPERSION")[0]

    e_bins = edisp["ENERG_LO"][1:]
    migra_bins = edisp["MIGRA_LO"][1:]

    fig = plt.gca()
    plt.pcolormesh(
        e_bins.to_value(u.TeV),
        migra_bins,
        edisp["MATRIX"].T[1:-1, 1:-1, 0].T,
        cmap="inferno",
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.colorbar(label="PDF Value")

    plt.xlabel(r"$E_\mathrm{True} / \mathrm{TeV}$")
    plt.ylabel(r"$E_\mathrm{Reco} / E_\mathrm{True}$")

    return fig


def plot_energy_resolution_from_file(filename, ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax

    bias_resolution = QTable.read(filename, hdu="ENERGY_BIAS_RESOLUTION")[1:-1]

    kwargs.setdefault("ls", "")

    # Plot function
    ax.errorbar(
        0.5 * (bias_resolution["true_energy_low"] + bias_resolution["true_energy_high"]).to_value(u.TeV),
        bias_resolution["resolution"],
        xerr=0.5 * (bias_resolution["true_energy_high"] - bias_resolution["true_energy_low"]).to_value(u.TeV),
        **kwargs,
    )
    ax.set_xscale("log")

    # Style settings
    ax.set_xlabel(r"$E_\mathrm{True} / \mathrm{TeV}$")
    ax.set_ylabel("Energy resolution")
    ax.grid(True, which="both")

    return ax


def plot_energy_bias_from_file(filename, ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax

    bias_resolution = QTable.read(filename, hdu="ENERGY_BIAS_RESOLUTION")[1:-1]

    kwargs.setdefault("ls", "")

    # Plot function
    ax.errorbar(
        0.5 * (bias_resolution["true_energy_low"] + bias_resolution["true_energy_high"]).to_value(u.TeV),
        bias_resolution["bias"],
        xerr=0.5 * (bias_resolution["true_energy_high"] - bias_resolution["true_energy_low"]).to_value(u.TeV),
        **kwargs,
    )
    ax.set_xscale("log")

    # Style settings
    ax.set_xlabel(r"$E_\mathrm{True}$ / $\mathrm{TeV}$")
    ax.set_ylabel("Energy bias")
    ax.grid(True, which="both")

    return ax

