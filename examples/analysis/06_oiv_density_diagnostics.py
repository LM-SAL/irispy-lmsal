"""
==================================
Measure densities with O IV ratios
==================================

In this example we measure the O IV 1399.78 / 1401.16 Å line ratio in a flare raster
and map it onto theoretical density curves with
`~irispy.utils.density.density_diagnostic`.

.. warning::

    This example requires ``fiasco`` and its CHIANTI database.
"""

import fiasco
import matplotlib.pyplot as plt
import numpy as np
import pooch

import astropy.units as u

from irispy.io import read_files
from irispy.utils.density import density_diagnostic
from irispy.utils.moments import average_window
from irispy.utils.spectrograph import subtract_background

###############################################################################
# The ratio of the O IV intercombination lines at 1399.78 and 1401.16 Å depends on the
# electron density :cite:p:`dudik2014,polito2016`. The O IV line at 1404.81 Å is blended
# with S IV 1404.85 Å :cite:p:`polito2016`, so we leave it out.
#
# We draw the curves at the three temperatures of Fig. 4 of :cite:t:`dudik2014`.
# `fiasco` computes the contribution function of every O IV transition at each
# temperature and density. A line is the sum of the transitions within 0.02 Å
# of its wavelength, as O IV 1401.149 and 1401.163 Å make up the 1401.157 Å line.

density = np.logspace(8, 13, 26) * u.cm**-3
log_temperatures = [4.80, 5.15, 5.50]
line_styles = [":", "-", "--"]
ion = fiasco.Ion("O IV", 10 ** np.array(log_temperatures) * u.K, ask_before=False)
contribution = ion.contribution_function(density, use_two_ion_model=False)
wavelengths = ion.transitions.wavelength[ion.transitions.is_bound_bound]


def line(wavelength):
    # The next O IV transitions lie 0.19 Å or more from these lines in CHIANTI 11.0.2
    return contribution[..., np.abs(wavelengths - wavelength) < 0.02 * u.AA].sum(axis=-1)


theoretical_ratio = line(1399.776 * u.AA) / line(1401.157 * u.AA)

###############################################################################
# We use the flare raster of :ref:`sphx_glr_generated_gallery_how_to_06_identify_lines.py`,
# and read only its Si IV 1394 window, with ``uncertainty=True`` so that the line
# intensities carry their uncertainties.

raster_filename = pooch.retrieve(
    "https://github.com/LM-SAL/irispy-data/releases/download/v1/iris_l2_20150130_055150_3893010094_cutout_raster.fits.gz",
    known_hash="603aa2a5dbe0cf9738e3628451dd24361da05b8eeecd42962ac1253db05888eb",
)
fuv2 = read_files(raster_filename, spectral_windows="Si IV 1394", uncertainty=True)["Si IV 1394"][0]

###############################################################################
# A constant fitted between the two lines, at 1400.2-1400.9 Å, removes the continuum.
# Both lines come from the same ion, so we sum each over the same Doppler band. It runs
# from -40 to +60 km/s to take in their red wings and stops short of S I 1401.51 Å
# :cite:p:`polito2016`.

fuv2 = subtract_background(fuv2, [1400.2, 1400.9] * u.AA, degree=0)
band = [-40, 60] * u.km / u.s
intensity_1399, intensity_1401 = [
    average_window(fuv2, band.to(u.AA, equivalencies=u.doppler_optical(rest)), method="sum")
    for rest in [1399.776, 1401.157] * u.AA
]

###############################################################################
# `~irispy.utils.density.density_diagnostic` maps the ratio of every pixel onto the
# curve at :math:`\log(T / \mathrm{K}) = 5.15`, the temperature of maximum abundance of
# O IV :cite:p:`polito2016,dere2023`, and gives the densities at the ratio plus and minus
# its uncertainty as bounds. We select pixels by the signal-to-noise ratio (SNR) of the
# ratio. These uncertainties are those of the samples: the error of the fitted background
# constant, shared by every sample of a band, is not included and would make the
# uncertainty of the ratio about 12% larger.

observed = density_diagnostic(intensity_1399, intensity_1401, density, theoretical_ratio[1])
snr = observed["ratio"] / observed["ratio_uncertainty"]
on_curve = np.isfinite(observed["density"])
above_limit = observed["ratio"] > theoretical_ratio[1].max()
for minimum_snr in [3, 10]:
    selected = snr > minimum_snr
    print(
        f"SNR > {minimum_snr}: {selected.sum()} pixels, {(selected & on_curve).sum()} on the curve, "
        f"{(selected & above_limit).sum()} above its high-density limit"
    )

###############################################################################
# We plot the pixels with SNR > 10 on the curves, with the uncertainty of their ratio
# and their density bounds.

good = (snr > 10) & on_curve
log_density = np.log10(observed["density"][good].to_value(u.cm**-3))
log_bounds = [np.log10(observed[key][good].to_value(u.cm**-3)) for key in ["density_lower", "density_upper"]]

fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
for line_style, log_temperature, ratio in zip(line_styles, log_temperatures, theoretical_ratio, strict=True):
    ax.plot(
        np.log10(density.to_value(u.cm**-3)),
        ratio.value,
        color="black",
        linestyle=line_style,
        linewidth=1.8,
        label=rf"$\log(T / \mathrm{{K}}) = {log_temperature:.2f}$",
    )
ax.errorbar(
    log_density,
    observed["ratio"][good].value,
    xerr=[log_density - log_bounds[0], log_bounds[1] - log_density],
    yerr=observed["ratio_uncertainty"][good].value,
    color="tab:red",
    marker="o",
    markersize=3,
    linestyle="none",
    elinewidth=0.5,
    label=f"{good.sum()} pixels with SNR > 10",
)
ax.set_title("O IV 1399.78 / 1401.16 Å")
ax.set_xlabel(r"$\log(n_e / \mathrm{cm^{-3}})$")
ax.set_ylabel("Line intensity ratio")
ax.legend(loc="upper left")

plt.show()

###############################################################################
# Nine in ten of the pixels with SNR > 10 lie between about :math:`10^{10.5}` and
# :math:`10^{11.3}` cm\ :sup:`-3`. Nearly a quarter of the pixels with SNR > 3, and
# none with SNR > 10, have a ratio above the high-density limit of the curve, about 0.42,
# where `~irispy.utils.density.density_diagnostic` gives NaN; all but 19 of them lie
# within twice their uncertainty of it.
