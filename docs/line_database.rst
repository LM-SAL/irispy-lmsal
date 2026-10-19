.. _irispy-line-database:

***************************
IRIS spectral line database
***************************

``irispy`` includes an offline catalog of the atomic lines in the three IRIS spectrograph passbands, taken from the NIST Atomic Spectra Database and CHIANTI, with predicted strengths in a few reference solar atmospheres.
One can query it with `irispy.utils.lines.get_lines`.

These are the five strongest predicted FUV2 lines in the quiet Sun with coronal abundances:

.. code-block:: python

    >>> import astropy.units as u
    >>> from irispy.utils.lines import get_lines

    >>> fuv2 = get_lines([1389, 1407] * u.angstrom, region="quiet_sun", abundance="coronal")
    >>> for line in fuv2[:5]:
    ...     print(f"{line['ion']:6} {line['wavelength'].value:8.2f} {line['intensity_quiet_sun_coronal']:5.2f}  {line['wavelength_source']}")
    Si IV   1393.76  1.00  observed
    Si IV   1402.77  0.50  observed
    O IV    1401.16  0.29  ritz
    O IV    1404.81  0.08  chianti
    O IV    1399.78  0.07  chianti

Wavelengths are in vacuum.
Above 2000 Å most line lists give air wavelengths, which for Mg II are 0.82 Å (88 km/s) shorter (`Morton 2000 <https://doi.org/10.1086/317349>`__).
To find the window that holds a line, compare its wavelength with ``raster[window][0].meta.spectral_range`` for each window.
A wavelength range can be given in any spectral unit, and both ends are included.
Without ``region`` and ``abundance``, lines come in wavelength order.
With them, lines are ranked strongest first within each passband (FUV1, FUV2, then NUV), and only the documented lines and the stage III and higher lines at least 2% of the passband's strongest are kept.
Pass ``minimum_strength=0`` to rank every line; lines without a prediction come last.
``main_only=True`` returns only the :ref:`documented IRIS lines <documented-lines>`, and ``categories`` filters by :ref:`category <line-categories>`.

The :ref:`gallery example <sphx_glr_generated_gallery_how_to_06_identify_lines.py>` marks candidate lines on an observation.

Passbands
---------

The catalog covers the passbands in Table 2 of :cite:t:`depontieu2014`: FUV1 1331.7--1358.4 Å, FUV2 1389.0--1407.0 Å, and NUV 2782.7--2835.1 Å.
FUV2 data extend down to about 1380.7 Å, but lines below 1389 Å are not in the catalog.

Columns
-------

* ``wavelength`` is NIST's observed or Ritz wavelength, whichever has the smaller uncertainty, as recorded in ``wavelength_source``.
  Both are kept in ``observed_wavelength`` and ``ritz_wavelength``.
  NIST notes that `Ritz wavelengths are usually more accurate in the vacuum ultraviolet <https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html>`__.
* Lines that are only in CHIANTI have ``wavelength_source="chianti"``.
  ``wavelength_is_theoretical`` marks the ones computed from theoretical energy levels, which can be less accurate: see `CHIANTI's wavelength conventions <https://db.chiantidatabase.org/o/o_4.html>`__.
  Documented lines in neither database have their published wavelength and ``wavelength_source="literature"``.
* ``lower`` and ``upper`` give each level's configuration, term, and J; ``lower_energy`` and ``upper_energy`` are in cm\ :sup:`-1`.
* ``nist_intensity`` is NIST's laboratory intensity, which depends on the light source and has no common scale.
* ``intensity_<region>_<abundance>`` are the eight :ref:`predicted strengths <line-strengths>`, NaN where CHIANTI has no prediction.
* ``strength_caveat`` is set for neutral and singly ionized ions, whose strengths are rough (see :ref:`line-limitations`).
* ``log_t_max`` is the log\ :sub:`10` temperature in K at which the ion's equilibrium fraction peaks.
* ``category`` is the ion's :ref:`category <line-categories>`.
* ``main`` marks the documented IRIS lines, and ``reference`` is their citation key.

The table's ``meta`` records the NIST and CHIANTI versions, the NIST queries, and the model settings used to build it.

.. _documented-lines:

Documented IRIS lines
---------------------

``main_only=True`` returns these published IRIS identifications.
Other lines in the catalog can matter too, such as the Si II and Fe II blends near Fe XXI in flares :cite:p:`young2015`.

======== ========== ===========================
Ion      Wavelength Reference
======== ========== ===========================
C II     1334.53    :cite:t:`depontieu2014`
C II     1335.66    :cite:t:`rathore2015`
C II     1335.71    :cite:t:`depontieu2014`
Fe XII   1349.40    :cite:t:`depontieu2014`
Cl I     1351.66    :cite:t:`itn38`
Fe XXI   1354.08    :cite:t:`depontieu2014`
C I      1354.29    :cite:t:`young2015`
Mn XVIII 1355.01    :cite:t:`depontieu2021`
O I      1355.60    :cite:t:`depontieu2014`
C I      1355.85    :cite:t:`lin2017`
Fe II    1392.82    :cite:t:`wulser2018`
Ni II    1393.33    :cite:t:`itn38`
Si IV    1393.76    :cite:t:`depontieu2014`
O IV     1399.78    :cite:t:`depontieu2014`
O IV     1401.16    :cite:t:`depontieu2014`
Si IV    1402.77    :cite:t:`depontieu2014`
O IV     1404.81    :cite:t:`polito2016`
S IV     1404.81    :cite:t:`polito2016`
S IV     1406.01    :cite:t:`polito2016`
Mg II    2791.60    :cite:t:`pereira2015`
Mg II    2796.35    :cite:t:`depontieu2014`
Mg II    2798.75    :cite:t:`pereira2015`
Mg II    2798.82    :cite:t:`pereira2015`
Ni I     2799.47    :cite:t:`wulser2018`
Mg II    2803.53    :cite:t:`depontieu2014`
======== ========== ===========================

O IV 1399.78, O IV 1404.81, and Mn XVIII are not in NIST and have CHIANTI's wavelengths; Fe II 1392.82 and Ni II 1393.33 are in neither and have their references' wavelengths.

.. _line-strengths:

Predicted strengths
-------------------

Each strength is the optically thin intensity

.. math::

   I = \frac{1}{4\pi}\int G(T, n_e)\,\mathrm{DEM}(T)\,\mathrm{d}T,

where :math:`\mathrm{DEM}(T)=n_e n_H\,\mathrm{d}h/\mathrm{d}T` and the contribution function :math:`G` comes from `fiasco <https://fiasco.readthedocs.io>`__ with CHIANTI 11.0.2 atomic data.
Strengths are divided by the strongest stage III or higher line in the same passband and model, so Si IV 1393.76 is 1 in FUV2.

* **DEM:** CHIANTI's quiet-Sun and active-region DEMs from `Vernazza and Reeves (1978) <https://doi.org/10.1086/190539>`__, and its flare DEM from the 1973 August 9 flare of `Dere and Cook (1979) <https://doi.org/10.1086/157013>`__.
* **Density:** a constant electron pressure of :math:`n_e T = 3\times10^{15}` K cm\ :sup:`-3`; ``flare_1e16`` is the flare DEM at :math:`10^{16}` K cm\ :sup:`-3`.
* **Ionization:** CHIANTI's equilibrium ionization fractions, which depend on temperature only.
  CHIANTI 11's density-dependent and charge-transfer models, which can change strong UV lines by factors of two to five :cite:p:`dufresne2024`, are not used.
* **Abundances:** ``photospheric`` is `Asplund, Amarsi, and Grevesse (2021) <https://doi.org/10.1051/0004-6361/202140445>`__; ``coronal`` is the same set with low first-ionization-potential elements such as Mg, Si, and Fe raised by 0.5 dex.

.. _line-limitations:

Limitations
-----------

The region names select fixed reference atmospheres; nothing is fitted to your observation.
Strengths rank integrated emission within one passband and model, and are not peak intensities or count rates.
Confirm any candidate against the observed spectrum and published identifications.

* **Neutral and singly ionized ions:** the optically thin model is far off for them.
  It puts Mg II k at 10\ :sup:`5` to 10\ :sup:`7` times the strongest NUV coronal line and Cr II above Si IV in active regions, so ranked queries keep them only when they are documented.
  Mg II h and k and C II need optically thick radiative transfer :cite:p:`leenaarts2013,rathore_carlsson2015`, and O I 1355.6 is set by recombination and charge exchange with hydrogen :cite:p:`lin2015`.
* **Opacity and equilibrium:** Si IV can become optically thick in flares :cite:p:`kerr2019`, and dynamic plasma can be out of ionization equilibrium.
* **Density:** at :math:`10^5` K the flare model has :math:`n_e=10^{11}` cm\ :sup:`-3`, while IRIS flare ribbons have shown about :math:`10^{13}` cm\ :sup:`-3`, which changes the O IV and S IV intensities and ratios :cite:p:`polito2016`.
* **DEMs:** they are 1970s averages derived with older coronal abundances, least certain at low temperatures, and the same DEM is used with both abundance sets.
  The two abundance sets show how sensitive a strength is to abundance, not how accurate it is.
* **Not modeled:** absorption, photoexcitation, sunspots, and molecular lines such as H\ :sub:`2` and CO.

.. _line-categories:

Categories
----------

Categories describe ions, not where or how a particular line forms.
The temperature boundaries follow :cite:t:`tian2017,depontieu2021` and use ``log_t_max``:

* ``flare``: :math:`T_\mathrm{max} \geq 10` MK, such as Fe XXI.
* ``coronal``: :math:`0.8 \leq T_\mathrm{max} < 10` MK, such as Fe XII :cite:p:`testa2016` and Mn XVIII.
* ``transition_region``: :math:`0.02 \leq T_\mathrm{max} < 0.8` MK, such as Si IV, O IV, and S IV.
* ``chromospheric``: Mg II, C II, O I, C I, and Cl I, following their IRIS formation studies :cite:p:`leenaarts2013,rathore_carlsson2015,lin2015,lin2017,itn38`.
* ``cool_metal``: other neutral and singly ionized metals from lithium to zinc, such as Si II, Fe I, and Fe II.
  These include photospheric absorption lines and chromospheric emission, such as the blends near Fe XXI in flares :cite:p:`young2015`.

Other lines, such as laboratory lines of heavy elements, have an empty category.
