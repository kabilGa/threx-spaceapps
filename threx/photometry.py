# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""How bright would the body be?  Reflected sunlight, plus a check on its own heat glow.

Plain words: the planet has no light of its own at these wavelengths. It is a mirror ball
lit by the Sun. Sunlight fades as 1/r^2 on the way out, the planet catches a part set by its
size R and reflects a fraction p (albedo), and the light fades again as 1/Delta^2 on the way
back. So brightness drops as 1/(r^2 Delta^2), roughly distance to the 4th power: twice as
far means 16 times fainter (3 magnitudes).
"""
import numpy as np

from .constants import AU_KM, R_EARTH_KM

# Apparent AB magnitude of the Sun at 1 AU in standard filters.
# MEASURED/tabulated: Willmer (2018), ApJS 236, 47, online table
# http://mips.as.arizona.edu/~cnaw/sun.html (column "apparent AB magnitude"),
# read on 2026-11-15. (pivot wavelength in micron, magnitude).
SUN_AB_1AU_TABLE = [
    (0.5511, -26.77, "Johnson V"),
    (0.7490, -27.05, "SDSS i"),
    (0.8947, -27.07, "SDSS z"),
    (1.2393, -27.03, "2MASS J"),
    (1.6495, -26.91, "2MASS H"),
    (2.1638, -26.49, "2MASS Ks"),
    (3.3897, -25.66, "WISE W1"),
    (3.5508, -25.56, "IRAC 3.6"),
    (3.8247, -25.40, "Bessell L'"),
    (4.4960, -25.06, "IRAC 4.5"),
    (4.6406, -25.00, "WISE W2"),
    (4.7347, -24.93, "Bessell M"),
]
SUN_V_AB_1AU = -26.77

# Albedo scenarios. ALL ASSUMED (nobody has seen this planet).
# - 0.2 to 0.6: the range asked for in our brief; Fortney et al. 2016 (arXiv:1604.07424)
#   argue a cold Planet Nine could even reach ~0.75 in visible light (Rayleigh scattering).
# - "methane-dark": Uranus and Neptune reflect only a few percent in the methane bands of the
#   near infrared. 0.03 is an ILLUSTRATIVE order of magnitude, not a fitted value.
ALBEDO_SCENARIOS = {
    "Bright (p = 0.6)": 0.6,
    "Middle (p = 0.4)": 0.4,
    "Dim (p = 0.2)": 0.2,
    "Methane-dark near-IR (p = 0.03, illustrative)": 0.03,
}
# Radius is ASSUMED from the mass (about 6 Earth masses): 2 to 3.7 Earth radii in our brief;
# Fortney et al. 2016 quote 3 to 6 Earth radii depending on the hydrogen/helium envelope.
RADIUS_RANGE_REARTH = (2.0, 3.7)


def sun_ab_mag(lam_um):
    """Apparent AB magnitude of the Sun at 1 AU at wavelength lam_um (micron).

    Interpolated (linear in log wavelength) between the Willmer (2018) filter values.
    Beyond 4.73 micron we extend with the Rayleigh-Jeans slope (flux ~ 1/lambda^2).
    """
    lam = np.asarray(lam_um, dtype=float)
    x = np.log10([row[0] for row in SUN_AB_1AU_TABLE])
    y = np.array([row[1] for row in SUN_AB_1AU_TABLE])
    m = np.interp(np.log10(lam), x, y)
    last = SUN_AB_1AU_TABLE[-1][0]
    return np.where(lam > last, y[-1] + 5.0 * np.log10(lam / last), m)


def sun_ab_mag_blackbody(lam_um, teff=5772.0, r_sun_km=695_700.0):
    """Cross-check only: the Sun as a perfect 5772 K blackbody (LABELLED APPROXIMATION)."""
    f_jy = _blackbody_flux_jy(teff, r_sun_km / AU_KM, np.asarray(lam_um, dtype=float))
    return -2.5 * np.log10(f_jy / 3631.0)


def reflected_mag(r_au, delta_au, radius_rearth, albedo, lam_um=None, m_sun=None):
    """AB magnitude of reflected sunlight at zero phase angle.

        m = m_sun(1 AU) - 2.5 log10( p * (R / 1 AU)^2 / (r^2 * Delta^2) )

    r_au    : Sun-planet distance [AU]      delta_au : Earth-planet distance [AU]
    radius_rearth : planet radius [Earth radii]      albedo : geometric albedo p
    lam_um  : wavelength [micron] used to look up the Sun's magnitude (default: V band).
    """
    if m_sun is None:
        m_sun = SUN_V_AB_1AU if lam_um is None else sun_ab_mag(lam_um)
    R_au = np.asarray(radius_rearth, dtype=float) * R_EARTH_KM / AU_KM
    ratio = albedo * R_au ** 2 / (np.asarray(r_au, dtype=float) ** 2 * np.asarray(delta_au, dtype=float) ** 2)
    return m_sun - 2.5 * np.log10(ratio)


def _blackbody_flux_jy(temp_k, angular_radius_rad, lam_um):
    """Flux density (Jy) of a blackbody disc of given angular radius: F = pi B_nu theta^2.
    Computed through logarithms so that very cold bodies do not underflow to zero."""
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    nu = c / (np.asarray(lam_um, dtype=float) * 1e-6)
    x = h * nu / (k * temp_k)
    # log10 of B_nu in W m^-2 Hz^-1 sr^-1 ; for large x, 1/(e^x - 1) -> e^-x
    log_planck = np.where(x > 30, -x / np.log(10), -np.log10(np.expm1(np.minimum(x, 30))))
    log_b = np.log10(2 * h * nu ** 3 / c ** 2) + log_planck
    log_f = log_b + np.log10(np.pi * angular_radius_rad ** 2) + 26.0     # 1 Jy = 1e-26 W m^-2 Hz^-1
    return 10.0 ** np.clip(log_f, -300, 300)


def thermal_mag(temp_k, radius_rearth, delta_au, lam_um):
    """AB magnitude of the planet's own blackbody glow at temperature temp_k.

    Caveat: real cold giant-planet atmospheres are NOT blackbodies. Fortney et al. 2016
    find 3-5 micron thermal fluxes many orders of magnitude above the blackbody value
    (and still very faint). We only use this as a 'is it negligible?' check.
    """
    theta = radius_rearth * R_EARTH_KM / (delta_au * AU_KM)
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    nu = c / (np.asarray(lam_um, dtype=float) * 1e-6)
    x = h * nu / (k * temp_k)
    log_planck = np.where(x > 30, -x / np.log(10), -np.log10(np.expm1(np.minimum(x, 30))))
    log_f_jy = np.log10(2 * h * nu ** 3 / c ** 2) + log_planck + np.log10(np.pi * theta ** 2) + 26.0
    return -2.5 * (log_f_jy - np.log10(3631.0))


def equilibrium_temperature(r_au, bond_albedo=0.3):
    """Temperature of a fast-rotating body heated only by the Sun: T = 278.6 K (1-A)^0.25 / sqrt(r)."""
    return 278.6 * (1.0 - bond_albedo) ** 0.25 / np.sqrt(r_au)
