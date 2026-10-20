.. _irispy-tutorial-fitting:

**********************
Fitting spectral lines
**********************

``irispy`` does not provide line fitting as astropy already handles this.
What ``irispy`` adds, in `irispy.utils.fitting`, are starting models for IRIS lines:

* `~irispy.utils.fitting.line_model`: a Gaussian for each documented line you name, on a constant background, with two for the emission peaks of Mg II k and h.
* `~irispy.utils.fitting.starts_from_model`: that model, or any model fitted to one spectrum, started for every spectrum from its own background and peaks.

Build a model with one value per parameter and let `~irispy.utils.fitting.starts_from_model` give each spectrum its own: a compound model built with array parameters fits every spectrum to NaN without an error (`astropy#20554 <https://github.com/astropy/astropy/issues/20554>`__).

Fitting a cube
==============

We use the Si IV 1394 window of the sample raster, which is the sample data of the other tutorials, and a small part of its field of view:

.. code-block:: python

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.coordinates import SpectralCoord
    >>> from astropy.modeling.fitting import TRFLSQFitter, parallel_fit_dask
    >>> import irispy.data.sample as sample_data  # doctest: +REMOTE_DATA
    >>> from irispy.io import read_files
    >>> from irispy.utils.fitting import line_model, starts_from_model

    >>> raster = read_files(sample_data.RASTER_FITS, spectral_windows="Si IV 1394", uncertainty=True)  # doctest: +REMOTE_DATA
    >>> cube = raster["Si IV 1394"][0]  # doctest: +REMOTE_DATA

`~irispy.utils.fitting.line_model` names the lines by their ion and wavelength in Å, as in the :doc:`line catalog </line_database>`.
The model has only the lines that you name, so the cube has to be cropped to them.
The Si IV 1394 window holds Fe II, Ni II and O IV lines as well as Si IV 139.376 nm; the catalog lists them.
Here we keep the wavelengths halfway to the Ni II and O IV lines on either side of Si IV:

.. code-block:: python

    >>> from irispy.utils.lines import get_lines
    >>> lines = get_lines([1390, 1400] * u.AA, main_only=True)
    >>> print(", ".join(lines["ion"]))
    Fe II, Ni II, Si IV, O IV
    >>> ni_ii, si_iv, o_iv = lines["wavelength"][1:]
    >>> window = [(ni_ii + si_iv) / 2, (si_iv + o_iv) / 2]
    >>> cube = cube.crop([SpectralCoord(window[0]), None], [SpectralCoord(window[1]), None])  # doctest: +REMOTE_DATA
    >>> cube = cube[:, 410:470]  # doctest: +REMOTE_DATA

`~irispy.utils.fitting.line_model` puts a Gaussian at the rest wavelength of Si IV, with a standard deviation of 0.05 Å, on a constant background.
`~irispy.utils.fitting.starts_from_model` starts it in every spectrum at that spectrum's peak above its background:

.. code-block:: python

    >>> model = starts_from_model(line_model(["Si IV 1394"]), cube)  # doctest: +REMOTE_DATA
    >>> model.param_names  # doctest: +REMOTE_DATA
    ('amplitude_0', 'amplitude_1', 'mean_1', 'stddev_1')
    >>> model.mean_1.shape  # doctest: +REMOTE_DATA
    (16, 60)

The background is component 0 and the line component 1.
Now we fit every spectrum:

.. code-block:: python

    >>> wavelength = cube.axis_world_coords("em.wl")[0].to(u.AA)  # doctest: +REMOTE_DATA
    >>> good = np.isfinite(cube.data) & ~cube.mask  # doctest: +REMOTE_DATA
    >>> fitter = TRFLSQFitter()
    >>> fitted = parallel_fit_dask(
    ...     model=model,
    ...     fitter=fitter,
    ...     data=np.where(good, cube.data, 0),
    ...     data_unit=cube.unit,
    ...     weights=np.where(good, 1 / cube.uncertainty.array, 0),
    ...     world=(wavelength,),
    ...     fitting_axes=2,
    ...     fit_info=True,
    ...     scheduler="single-threaded",
    ... )  # doctest: +REMOTE_DATA
    >>> covariance = fitter.fit_info.get_property_as_array("param_cov")  # doctest: +REMOTE_DATA
    >>> errors = np.sqrt(np.diagonal(covariance, axis1=-2, axis2=-1))  # doctest: +REMOTE_DATA
    >>> fitted.mean_1.quantity.shape, errors.shape  # doctest: +REMOTE_DATA
    ((16, 60), (16, 60, 4))

``errors[..., i]`` is the 1σ uncertainty of the parameter ``fitted.param_names[i]``.
The other properties of each fit, such as ``status``, ``success`` and ``nfev``, come out of ``fitter.fit_info`` in the same way.
The default ``scheduler`` uses several processes; ``"single-threaded"`` keeps this example simple.
A ``dask.distributed.Client`` can be passed instead.
Any dask scheduler works. Pickling a cube is not a supported way to save or share results; re-run the script that made them.

Maps from the fit
=================

`~irispy.utils.fitting.maps_from_fit` turns the fitted model into a `~irispy.spectrograph.RasterCollection` of maps with the spatial coordinates of the cube.
When the fitter provides covariance, parameter maps carry the uncertainties from the fit, and derived maps propagate them:

.. code-block:: python

    >>> from irispy.utils.fitting import maps_from_fit
    >>> maps = maps_from_fit(fitted, cube, fitter=fitter)  # doctest: +REMOTE_DATA
    >>> list(maps.keys())  # doctest: +REMOTE_DATA
    ['amplitude_0', 'amplitude_1', 'mean_1', 'stddev_1', 'velocity_1', 'fwhm_1', 'fwhm_velocity_1', 'integrated_intensity_1', 'quality', 'residual']
    >>> maps["velocity_1"].unit, maps["integrated_intensity_1"].unit  # doctest: +REMOTE_DATA
    (Unit("km / s"), Unit("Angstrom DN_IRIS_FUV"))

The velocities are against the rest wavelength of the one documented line in the cube, here Si IV, unless ``rest_wavelength`` is given.
``"quality"`` holds a `~irispy.utils.fitting.FitQualityFlag` for each spectrum: whether its fit failed or did not converge, stopped at a bound, or used a spectrum with masked samples.
``"residual"`` is the cube minus the fitted model, to check the fits. It retains the cube's uncertainty and coordinates.

Judging the fits
----------------

`~irispy.utils.fitting.FitQualityFlag.OK` only says that the fitter converged within the bounds of the parameters.
A fit to noise converges just as well, so it does not say that there is a line.
Two checks, both from quantities above, tell more.
The signal-to-noise ratio of the line is the fitted amplitude over its 1σ error, and the reduced chi-square is the sum of the squared residuals over the uncertainties of the unmasked samples, divided by their number minus the number of parameters:

.. code-block:: python

    >>> from irispy.utils.fitting import FitQualityFlag
    >>> with np.errstate(divide="ignore", invalid="ignore"):  # doctest: +REMOTE_DATA
    ...     snr = fitted.amplitude_1.value / errors[..., 1]
    >>> chi2 = np.sum(np.where(good, maps["residual"].data / cube.uncertainty.array, 0) ** 2, axis=-1)  # doctest: +REMOTE_DATA
    >>> reduced_chi2 = chi2 / (good.sum(axis=-1) - len(fitted.param_names))  # doctest: +REMOTE_DATA
    >>> converged = maps["quality"].data == FitQualityFlag.OK  # doctest: +REMOTE_DATA
    >>> f"{converged.mean():.0%}", int((converged & (snr > 3)).sum())  # doctest: +REMOTE_DATA
    ('80%', 73)

Here the fitter converged for 80% of the 960 spectra, but only 73 of those have a signal-to-noise ratio above 3.
There is no universal threshold: choose one for your science case, and keep the fits above it.
A reduced chi-square near 1 means that the model describes the spectrum within its uncertainties, which a fit to pure noise does as well, so it is a check of the model and of the uncertainties and not of the presence of the line.
Values well above 1 point to a profile that the model does not describe, such as a blend or a second component, or to uncertainties that are too small.

Fitting the mean spectrum of a region
=====================================

A line too faint to fit in each spectrum can be fitted in the mean spectrum of a region, as ``iris_sum_spec`` and ``dgf_1lp`` do in IDL, and that fit can start the fits of every spectrum.
Here the region is the whole raster.
Cosmic rays are not removed from Level 2 data and a mean takes them in, so we remove them first with `~irispy.utils.cosmic_rays.remove_cosmic_rays`, which needs the ``cosmic-rays`` extra (see :ref:`irispy-install-extras`).
The mean leaves out masked and missing samples, and its uncertainty at each wavelength is :math:`\sqrt{\sum \sigma^2} / n` for the :math:`n` samples averaged:

.. code-block:: python

    >>> from irispy.utils.cosmic_rays import remove_cosmic_rays
    >>> region = raster["Si IV 1394"][0].crop([SpectralCoord(window[0]), None], [SpectralCoord(window[1]), None])  # doctest: +REMOTE_DATA
    >>> region = remove_cosmic_rays(region)  # doctest: +REMOTE_DATA
    >>> usable = np.isfinite(region.data) & ~region.mask  # doctest: +REMOTE_DATA
    >>> spectra = np.ma.masked_array(region.data, ~usable)  # doctest: +REMOTE_DATA
    >>> sigma = np.ma.masked_array(region.uncertainty.array, ~usable)  # doctest: +REMOTE_DATA
    >>> spectrum = spectra.mean(axis=(0, 1))  # doctest: +REMOTE_DATA
    >>> spectrum_error = np.sqrt(np.sum(sigma**2, axis=(0, 1))) / spectra.count(axis=(0, 1))  # doctest: +REMOTE_DATA

``region.rebin`` with ``propagate_uncertainties=True`` also averages with uncertainties, but its cost grows with the square of the number of pixels in a bin, so it gets slow for large regions.

The model, built with astropy, is a narrow core and a broad component on a constant background.
The fit may return them in either order, so the narrower one is the core:

.. code-block:: python

    >>> from astropy.modeling import models
    >>> core_and_broad = (
    ...     models.Const1D(0.3 * region.unit)
    ...     + models.Gaussian1D(1 * region.unit, si_iv, 0.03 * u.AA)
    ...     + models.Gaussian1D(1 * region.unit, si_iv, 0.1 * u.AA)
    ... )  # doctest: +REMOTE_DATA
    >>> mean_fit = TRFLSQFitter()(core_and_broad, wavelength, spectrum * region.unit, weights=1 / spectrum_error)  # doctest: +REMOTE_DATA
    >>> core, broad = sorted([1, 2], key=lambda i: mean_fit[i].stddev)  # doctest: +REMOTE_DATA
    >>> print(f"{mean_fit[core].stddev.quantity:.2f}, {mean_fit[broad].stddev.quantity:.2f}")  # doctest: +REMOTE_DATA
    0.03 Angstrom, 0.10 Angstrom

We remove the cosmic rays from the cube as well, and `~irispy.utils.fitting.starts_from_model` starts each of its spectra from this fit.
These spectra are too faint to fit all seven parameters, so we fix the widths at those of the mean spectrum and give both components one center.
Those widths also hold the spread of Doppler shifts over the region, so they describe the region and not each spectrum:

.. code-block:: python

    >>> mean_fit.stddev_1.fixed = mean_fit.stddev_2.fixed = True  # doctest: +REMOTE_DATA
    >>> mean_fit.mean_2.tied = lambda model: model.mean_1  # doctest: +REMOTE_DATA
    >>> clean_cube = remove_cosmic_rays(cube)  # doctest: +REMOTE_DATA
    >>> clean_good = np.isfinite(clean_cube.data) & ~clean_cube.mask  # doctest: +REMOTE_DATA
    >>> cube_fit = parallel_fit_dask(
    ...     model=starts_from_model(mean_fit, clean_cube),
    ...     fitter=TRFLSQFitter(),
    ...     data=np.where(clean_good, clean_cube.data, 0),
    ...     data_unit=clean_cube.unit,
    ...     weights=np.where(clean_good, 1 / clean_cube.uncertainty.array, 0),
    ...     world=(wavelength,),
    ...     fitting_axes=2,
    ...     scheduler="single-threaded",
    ... )  # doctest: +REMOTE_DATA
    >>> cube_fit[core].amplitude.shape  # doctest: +REMOTE_DATA
    (16, 60)

Non-thermal velocities
======================

A Gaussian fitted width holds the instrumental and thermal broadening of the line as well.
The quadrature subtraction assumes Gaussian profiles; do not pass the width of a Lorentzian fit.

`~irispy.utils.fitting.non_thermal_velocity` removes both from a map of the full width at half maximum, using the spectral resolution of the passband and, for the thermal width, the ion's mass and a temperature, here the peak of Si IV in CHIANTI's ionization equilibrium :cite:p:`dere2023`:

.. code-block:: python

    >>> from irispy.utils.fitting import non_thermal_velocity
    >>> non_thermal = non_thermal_velocity(maps["fwhm_1"], si_iv, ion="Si IV", temperature=10**4.9 * u.K)  # doctest: +REMOTE_DATA
    >>> non_thermal.unit  # doctest: +REMOTE_DATA
    Unit("km / s")

The velocity is masked where the width is masked or missing, or not above the other two;
``non_thermal.mask & ~maps["fwhm_1"].mask & np.isfinite(maps["fwhm_1"].data)`` picks out the pixels where it is too narrow.

Things to know
==============

* **Pass the wavelengths in Å**
  The starting models are in Å, and astropy converts a model's parameters to the unit of the wavelengths it is given but not their bounds.
  `~astropy.modeling.fitting.TRFLSQFitter` also does not rescale the parameters, so with wavelengths in meters, the unit of the cube's WCS, it can stop before the line centers and widths converge.
  That is why the wavelengths are passed as ``world`` rather than the cube itself.
* **Masked and missing samples**
  The fitters do not accept NaN, so replace bad samples with any number and give them zero weight, as above.
  Fits that raise an error give NaN parameters and a zero covariance, and spectra left entirely NaN make ``fit_info.get_property_as_array`` fail; `~irispy.utils.fitting.maps_from_fit` copes with both.
  ``diagnostics="error"`` with a ``diagnostics_path`` writes each error to a folder.
* **Uncertainties**
  With weights of 1/σ, ``param_cov`` is the covariance of the parameters given those uncertainties, and the errors are only as good as the uncertainties of the cube.
  Without weights, astropy scales it by the variance of the residuals, as `scipy.optimize.curve_fit` does with ``absolute_sigma=False``.
* **Bounds**
  `~astropy.modeling.fitting.TRFLSQFitter` keeps the parameters within their bounds while it fits;
  `~astropy.modeling.fitting.LMLSQFitter`, which is faster, clips them to the bounds at each step.
* **Faint lines**
  When a line is not much brighter than the noise, its fit is poorly constrained whatever its start.
  Average neighboring spectra first, keeping the uncertainties with ``cube.rebin((2, 2, 1), propagate_uncertainties=True)``.
* **Two Gaussians can swap**
  In `~irispy.utils.fitting.line_model`, the k2v Gaussian starts at the blue peak and the k2r one at the red peak, but both may move anywhere within 40 km/s of the rest wavelength.
  A fit can end with them the other way round, so sort the components by velocity before combining their maps.
* **Several components**
  `~irispy.utils.fitting.line_model` gives each line one Gaussian because the layout of several depends on what is observed:
  a narrow and a broad Gaussian with one centroid in active-region loops :cite:p:`dudik2017`, a static and a redshifted component in flare ribbons :cite:p:`yu2020`, and components tens of km/s from the line center in UV bursts :cite:p:`peter2014,young2018`.
  Build such a model with astropy, as for the mean spectrum above, and start it with `~irispy.utils.fitting.starts_from_model`.

The gallery examples :ref:`sphx_glr_generated_gallery_analysis_01_spectral_fitting.py` and :ref:`sphx_glr_generated_gallery_analysis_07_mg_ii_two_gaussian_fitting.py` plot such maps.
