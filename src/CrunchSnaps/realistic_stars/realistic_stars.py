"""Routines for rendering stars with realistic Hubble PSFs

Most credit goes to Pelupessy & Riedier, authors of the AMUSE fresco package. 

Ported here to avoid an AMUSE dependency

See https://github.com/rieder/amuse-fresco for the original, much more feature-rich code.
"""

from astropy.io import fits
from astropy import units as u, constants as c
import os
import numpy as np
from meshoid import Meshoid
from meshoid.radiation import dust_ext_opacity, radtransfer
from scipy.signal import convolve
from matplotlib import pyplot as plt
from starforge_tools import star_gas_columns
from starforge_tools.special_functions import planck_wavelength_integral


FILTER_WAVELENGTHS_NM = {"u": 365, "b": 445, "v": 551, "r": 658, "i": 806}
FILTER_WIDTHS_NM = {"u": 66, "b": 94, "v": 88, "r": 138, "i": 149}


def get_psf(instrument="WFPC_II_WFC3", sample_fac: int = 1) -> dict:
    """
    Returns a dict containing the PSF images for the specified instrument.
    Keys are "<filter><number>" where the filters are ubvri and numbers run
    from 0 to 3.

    Parameters
    ----------
    instrument: str, optional
        Hubble instrument; currently on WFPC_II_WFC3 is available.

    Returns
    -------
    psf: dict
        Dictionary containing PSFs
    """
    psf = dict()
    this_dir = os.path.dirname(os.path.abspath(__file__))
    datadir = this_dir + "/psf/" + instrument + "/"

    NUM_FILTER_PSFS = 4

    for band in FILTER_WAVELENGTHS_NM.keys():
        psf[band] = 0
        for i in range(NUM_FILTER_PSFS):
            psf_filename = datadir + band + "%2.2i.fits" % i
            f = fits.open(psf_filename)
            psf[band + str(i)] = np.array(f[0].data[::sample_fac, ::sample_fac])
            psf[band] += psf[band + str(i)] / NUM_FILTER_PSFS
    return psf


def lum_radius_to_Teff(lum, radius):
    return 5770 * np.sqrt(np.sqrt((lum / radius**2)))


def make_stars_image_fullRT(snapdata, lum_max_solar=1e3, IMG_RES: int = 2048, IMG_SIZE: float = 10.0, verbose=False):
    PIXEL_SIZE = IMG_SIZE / IMG_RES
    xstar = snapdata["PartType5/Coordinates"]
    Lstar = snapdata["PartType5/StarLuminosity_Solar"]
    mstar = snapdata["PartType5/BH_Mass"]
    num_stars = mstar.shape[0]
    boxsize = snapdata["BoxSize"]
    Rstar = snapdata["PartType5/ProtoStellarRadius_inSolar"]
    xgas = snapdata["PartType0/Coordinates"]
    mgas = snapdata["PartType0/Masses"]
    hgas = snapdata["PartType0/SmoothingLength"]

    Teff = lum_radius_to_Teff(Lstar, Rstar)
    Lband = get_stellar_lum_in_bands(Lstar.clip(0, lum_max_solar), Teff)

    num_wavelengths = len(FILTER_WAVELENGTHS_NM)

    # dust opacity in cgs converted to solar - evaluated at 555nm
    wavelengths_um = np.array([l / 1e3 for l in FILTER_WAVELENGTHS_NM.values()])
    kappa_dust_codeunits = dust_ext_opacity(wavelengths_um).to(u.pc**2 / c.M_sun).value
    kappa_gas = np.array(len(mgas) * [kappa_dust_codeunits])
    j_gas = np.zeros_like(kappa_gas)  # assume dust does not emit

    # have to get the star properties now
    xstar = xstar  # (load_from_snapshot("Coordinates", 4, ".", 600) - center) @ coordinate_basis
    h_star = np.repeat(PIXEL_SIZE, num_stars)
    j_star = np.ones((num_stars, num_wavelengths))
    for i, band in enumerate(FILTER_WAVELENGTHS_NM.keys()):
        j_star[:, i] = Lband[band] / mstar
    kappa_stars = np.zeros_like(j_star)

    # now combine all emissivities, opacities, masses, kernel lengths
    j_all = np.atleast_2d(np.concatenate([j_gas, j_star]))  # 2D because this has shape (num_particles, num_bands)
    kappa_all = np.atleast_2d(np.concatenate([kappa_gas, kappa_stars]))  # ditto
    kappa_all = kappa_all.clip(1e-100)  # avoid floating-point errors
    h_all = np.concatenate([hgas, h_star])
    m_all = np.concatenate([mgas, mstar])
    x_all = np.concatenate([xgas, xstar], axis=0)
    if verbose:
        print("raytracing...")
    I = radtransfer(j_all, m_all, kappa_all, x_all, h_all, IMG_RES, IMG_SIZE, center=np.zeros(3) + 0.5 * boxsize)
    if verbose:
        print("done!")
    image = {}
    psfs = get_psf()

    for i, band in enumerate(FILTER_WAVELENGTHS_NM.keys()):
        image[band] = convolve(I[:, ::-1, i].T, psfs[band], mode="same")

    image_rgb = np.empty((IMG_RES, IMG_RES, 3))
    image_rgb[:, :, 0] = image["r"]
    image_rgb[:, :, 1] = image["v"]
    image_rgb[:, :, 2] = image["b"]
    image_rgb /= image_rgb.mean() * 10
    image_rgb = image_rgb.clip(0, 1)
    plt.imsave("light_fullRT.png", image_rgb)
    return image_rgb


def make_stars_image_starsonly(snapdata, lum_max_solar=1e3, IMG_RES: int = 2048, IMG_SIZE: float = 10.0, verbose=False):
    PIXEL_SIZE = IMG_SIZE / IMG_RES
    xstar = snapdata["PartType5/Coordinates"]
    Lstar = snapdata["PartType5/StarLuminosity_Solar"]
    mstar = snapdata["PartType5/BH_Mass"]
    num_stars = mstar.shape[0]
    boxsize = snapdata["BoxSize"]
    Rstar = snapdata["PartType5/ProtoStellarRadius_inSolar"]
    xgas = snapdata["PartType0/Coordinates"]
    mgas = snapdata["PartType0/Masses"]
    hgas = snapdata["PartType0/SmoothingLength"]

    Teff = lum_radius_to_Teff(Lstar, Rstar)

    # dust opacity in cgs converted to solar - evaluated at 555nm
    wavelengths_um = np.array([l / 1e3 for l in FILTER_WAVELENGTHS_NM.values()])
    kappa_dust_codeunits = dust_ext_opacity(wavelengths_um).to(u.pc**2 / c.M_sun).value

    star_columns = star_gas_columns(xstar, xgas, mgas, hgas)

    tau_dust = kappa_dust_codeunits * star_columns[:, None]

    Lstar_in_bands = get_stellar_lum_in_bands(Lstar.clip(0, lum_max_solar), Teff)
    attenuation = np.exp(-tau_dust)

    M = Meshoid(xstar, kernel_radius=np.repeat(PIXEL_SIZE, num_stars))  # , n_jobs=1)

    image = {}
    psfs = get_psf()

    for i, band in enumerate(FILTER_WAVELENGTHS_NM.keys()):
        I = M.SurfaceDensity(
            Lstar_in_bands[band] * attenuation[:, i],
            center=np.zeros(3) + boxsize / 2,
            conservative=False,
            size=IMG_SIZE,
            res=IMG_RES,
        )
        image[band] = convolve(I[:, ::-1].T, psfs[band], mode="same")

    image_rgb = np.empty((IMG_RES, IMG_RES, 3))
    image_rgb[:, :, 0] = image["r"]
    image_rgb[:, :, 1] = image["v"]
    image_rgb[:, :, 2] = image["b"]
    image_rgb /= image_rgb.mean() * 10
    image_rgb = image_rgb.clip(0, 1)
    plt.imsave("light_stars.png", image_rgb)
    return image_rgb


def get_stellar_lum_in_bands(Lstar, Teff):
    lum_band = {}
    for band, wavelength in FILTER_WAVELENGTHS_NM.items():
        bandwidth = FILTER_WIDTHS_NM[band]
        wavelength_min = max(wavelength - 0.5 * bandwidth, 0)
        wavelength_max = wavelength + 0.5 * bandwidth
        frac = planck_wavelength_integral(wavelength_min / 1e3, wavelength_max / 1e3, Teff)
        lum_band[band] = frac * Lstar

    return lum_band
