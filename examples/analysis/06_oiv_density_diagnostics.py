"""
==============================
O IV Density-Diagnostic Curves
==============================

The goal of this example is to show how to compute and plot density-sensitive
O IV line ratios and map an observed ratio onto them with
`~irispy.utils.density.density_diagnostic`.
These ratios are sensitive to the electron density in the solar atmosphere, and
they are often used to diagnose conditions in the solar transition region.

.. warning::

    This example requires ``fiasco`` and its CHIANTI database.
"""

import fiasco
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u

from irispy.utils.density import density_diagnostic

###############################################################################
# We will reproduce aspects of the top row of Fig. 4 from
# :cite:t:`dudik2014`, which shows the O IV line
# ratios as a function of electron density for three different Maxwellian temperatures.
#
# `fiasco` computes the contribution function of every O IV transition at each
# temperature and density. A line is the sum of the transitions within 0.02 Å
# of its wavelength, as O IV 1401.149 and 1401.163 Å make up the 1401.157 Å line.

density = np.logspace(9, 12, 20) * u.cm**-3
temperature_samples = 10 ** np.array([4.80, 5.15, 5.50]) * u.K
line_styles = [":", "-", "--"]
temperature_labels = ["4.80", "5.15", "5.50"]
ion = fiasco.Ion("O IV", temperature_samples, ask_before=False)
contribution = ion.contribution_function(density, use_two_ion_model=False)
wavelengths = ion.transitions.wavelength[ion.transitions.is_bound_bound]


def line(wavelength):
    # The next O IV transitions lie 0.19 Å or more from these lines in CHIANTI 11.0.2
    return contribution[..., np.abs(wavelengths - wavelength) < 0.02 * u.AA].sum(axis=-1)


ratio_definitions = [
    ("O IV 1401.16 / 1404.78 Å", line(1401.157 * u.AA) / line(1404.806 * u.AA)),
    ("O IV 1404.78 / 1399.77 Å", line(1404.806 * u.AA) / line(1399.780 * u.AA)),
]

###############################################################################
# The exact curve values will differ somewhat from the paper because this example
# uses the current CHIANTI database through `fiasco` rather than the CHIANTI
# version used by :cite:t:`dudik2014`.
#
# `~irispy.utils.density.density_diagnostic` maps measured line intensities
# onto one of these curves. As a synthetic observation, we take the ratio of the
# :math:`\log(T / \mathrm{K}) = 5.15` curve at its middle density.

fig, axes = plt.subplots(
    ncols=2,
    figsize=(11, 4.2),
    sharex=True,
)

for ax, (title, ratios) in zip(axes, ratio_definitions, strict=True):
    for line_style, label, ratio in zip(line_styles, temperature_labels, ratios, strict=True):
        ax.plot(
            np.log10(density.to_value("cm-3")),
            ratio.value,
            color="black",
            linestyle=line_style,
            linewidth=1.8,
            label=rf"$\log(T / \mathrm{{K}}) = {label}$",
        )
    observed = density_diagnostic(ratios[1, density.size // 2] * 100 * u.ct, 100 * u.ct, density, ratios[1])
    # Plotted after the curves so that it comes last in the legend.
    ax.plot(
        np.log10(observed["density"].to_value("cm-3")),
        observed["ratio"].value,
        color="tab:red",
        marker="o",
        linestyle="none",
        label=r"synthetic observation ($\log(T / \mathrm{K}) = 5.15$)",
    )
    ax.set_title(title)
    ax.set_xlabel(r"$\log(n_e / \mathrm{cm^{-3}})$")
    ax.set_ylabel("Line intensity ratio")

axes[0].legend(loc="upper left")

plt.show()
