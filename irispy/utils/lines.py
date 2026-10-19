"""
Query the packaged IRIS spectral line database.
"""

from functools import cache

import numpy as np

import astropy.units as u
from astropy.table import QTable

from irispy.data import ROOTDIR

__all__ = ["get_lines"]

_REGIONS = ("quiet_sun", "active_region", "flare", "flare_1e16")
_ABUNDANCES = ("coronal", "photospheric")
_CATEGORIES = ("flare", "coronal", "transition_region", "chromospheric", "cool_metal")
_PASSBANDS = ("FUV1", "FUV2", "NUV")
# The share of a passband's strongest predicted line below which a ranked query leaves a line out
_MINIMUM_STRENGTH = 0.02


@cache
def _load_lines():
    # Keep empty strings as strings rather than masked values.
    return QTable.read(ROOTDIR / "iris_lines.ecsv.gz", format="ascii.ecsv", fill_values=[])


def get_lines(
    wavelength_range=None,
    *,
    region=None,
    abundance=None,
    minimum_strength=None,
    main_only=False,
    categories=None,
):
    """
    Query the packaged IRIS spectral line catalog.

    Parameters
    ----------
    wavelength_range : `~astropy.units.Quantity`, optional
        Two endpoints in any spectral unit, in either order; both are included.
        Wavelengths are in vacuum.
    region : `str`, optional
        Rank by predicted strength in the ``'quiet_sun'``, ``'active_region'``,
        ``'flare'`` or ``'flare_1e16'`` reference model, strongest first within each
        passband. ``'flare_1e16'`` is the flare model at an electron pressure of
        1e16 K cm-3 instead of 3e15. Requires ``abundance``. Without it, lines are
        in wavelength order.
    abundance : `str`, optional
        ``'coronal'`` or ``'photospheric'``. Requires ``region``.
    minimum_strength : `float`, optional
        With ``region``, keep the documented lines and the stage III and higher
        lines at least this fraction of the passband's strongest. Defaults to 0.02.
        Pass 0 for every line, with unpredicted lines last.
    main_only : `bool`, optional
        Only the documented IRIS lines; ``reference`` holds their citation keys.
    categories : `str` or iterable of `str`, optional
        Only ions in these categories: ``'flare'``, ``'coronal'``,
        ``'transition_region'``, ``'chromospheric'`` or ``'cool_metal'``.

    Returns
    -------
    `~astropy.table.QTable`
        A copy of the selected lines. Strengths (``intensity_<region>_<abundance>``)
        are relative to the strongest stage III or higher line in each passband.

    Notes
    -----
    Strengths are optically thin estimates from fixed reference atmospheres; see
    :doc:`/line_database` for what they assume.
    """
    if region is not None and region not in _REGIONS:
        msg = f"region must be one of {_REGIONS}."
        raise ValueError(msg)
    if abundance is not None and abundance not in _ABUNDANCES:
        msg = f"abundance must be one of {_ABUNDANCES}."
        raise ValueError(msg)
    if (region is None) != (abundance is None):
        msg = "Supply region and abundance together to rank lines."
        raise ValueError(msg)
    if minimum_strength is not None and region is None:
        msg = "minimum_strength needs region and abundance."
        raise ValueError(msg)
    if region is not None and minimum_strength is None:
        minimum_strength = _MINIMUM_STRENGTH
    table = _load_lines()
    selected = np.ones(len(table), dtype=bool)
    if wavelength_range is not None:
        endpoints = u.Quantity(wavelength_range)
        if endpoints.shape != (2,) or not np.all(np.isfinite(endpoints.value)) or np.any(endpoints.value <= 0):
            msg = "wavelength_range must contain two finite, positive spectral endpoints."
            raise ValueError(msg)
        # Widen by a few ulps so endpoints survive spectral-unit round trips.
        low, high = np.sort(endpoints.to_value(u.angstrom, equivalencies=u.spectral())) * (
            1 + 4 * np.finfo(float).eps * np.array([-1, 1])
        )
        wavelength = table["wavelength"].to_value(u.angstrom)
        selected &= (wavelength >= low) & (wavelength <= high)
    if main_only:
        selected &= table["main"]
    if categories is not None:
        categories = [categories] if isinstance(categories, str) else list(categories)
        if any(category not in _CATEGORIES for category in categories):
            msg = f"categories must be drawn from {_CATEGORIES}."
            raise ValueError(msg)
        selected &= np.isin(table["category"], categories)

    intensity_column = f"intensity_{region}_{abundance}" if region is not None else None
    if intensity_column and minimum_strength > 0:
        # The cool ions' strengths are far outside the model's range (Mg II k at 1e5 times the strongest NUV
        # coronal line), so only the documented ones make the short list
        strength = table[intensity_column]
        hot = (table["ion_stage"] >= 3) & np.isfinite(strength) & (strength >= minimum_strength)
        selected &= table["main"] | hot

    # Isolate cached metadata.
    result = table[selected].copy()
    if intensity_column:
        intensity = result[intensity_column]
        passband = np.array([_PASSBANDS.index(band) for band in result["passband"]])
        order = np.lexsort((result["wavelength"].value, -intensity, passband))  # NaN sorts last
        return result[order]
    return result
