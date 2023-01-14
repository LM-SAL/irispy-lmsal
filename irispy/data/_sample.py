from pathlib import Path
from urllib.parse import urljoin

from sunpy import log
from sunpy.data._sample import _download_sample_data
from sunpy.util.config import get_and_create_sample_dir
from sunpy.util.parfive_helpers import Downloader

_BASE_URLS = (
    "https://github.com/sunpy/data/raw/main/irispy-lmsal/",
    "https://github.com/sunpy/sample-data/raw/master/irispy-lmsal/",
    "http://data.sunpy.org/irispy-lmsal/",
)
_SAMPLE_DATA = {
    "AIA_1700": "aia_20140919_060030_1700_image_lev1.fits",
    "RASTER": "iris_l2_20211001_060925_3683602040_raster.tar.gz",
    "SJI_1330": "iris_l2_20211001_060925_3683602040_SJI_1330_t000.fits.gz",
    "SJI_1400": "iris_l2_20211001_060925_3683602040_SJI_1400_t000.fits.gz",
    "SJI_2796": "iris_l2_20211001_060925_3683602040_SJI_2796_t000.fits.gz",
    "SJI_2832": "iris_l2_20211001_060925_3683602040_SJI_2832_t000.fits.gz",
}
_SAMPLE_FILES = {v: k for k, v in _SAMPLE_DATA.items()}  # NOQA


def _download_sample_data(base_url, sample_files, overwrite):
    """
    Downloads a list of files.

    Parameters
    ----------
    base_url : str
        Base URL for each file.
    sample_files : list of tuples
        List of tuples that are (URL_NAME, SAVE_NAME).
    overwrite : bool
        Will overwrite a file on disk if True.

    Returns
    -------
    `parfive.Results`
        Download results. Will behave like a list of files.
    """
    dl = Downloader(overwrite=overwrite, progress=True, headers={"Accept-Encoding": "identity"})
    for url_file_name, fname in sample_files:
        url = urljoin(base_url, url_file_name)
        dl.enqueue_file(url, filename=fname)
    results = dl.download()
    return results


def _retry_sample_data(results, new_url_base):
    dl = Downloader(overwrite=True, progress=True, headers={"Accept-Encoding": "identity"})
    for err in results.errors:
        file_name = err.url.split("/")[-1]
        log.debug(f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}")
        new_url = urljoin(new_url_base, file_name)
        log.debug(f"Attempting redownload of {_SAMPLE_FILES[file_name]} using {new_url}")
        dl.enqueue_file(new_url, filename=err.filepath_partial)
    extra_results = dl.download()
    new_results = results + extra_results
    new_results._errors = extra_results._errors
    return new_results


def _handle_final_errors(results):
    for err in results.errors:
        file_name = err.url.split("/")[-1]
        log.debug(f"Failed to download {_SAMPLE_FILES[file_name]} from {err.url}: {err.exception}")
        log.error(f"Failed to download {_SAMPLE_FILES[file_name]} from all mirrors, the file will not be available.")
        log.error(err)


def download_sample_data(overwrite=False):
    """
    Download all sample data at once. This will overwrite any existing files.

    Parameters
    ----------
    overwrite : `bool`
        Overwrite existing sample data.
    """
    sampledata_dir = Path(get_and_create_sample_dir()).parent / Path("irispy")
    already_downloaded = []
    to_download = []
    for url_file_name in _SAMPLE_FILES.keys():
        fname = sampledata_dir / url_file_name
        if fname.exists() and not overwrite:
            already_downloaded.append(fname)
        else:
            to_download.append((url_file_name, fname))
    if to_download:
        results = _download_sample_data(_BASE_URLS[0], to_download, overwrite=overwrite)
    else:
        return already_downloaded
    if results.errors:
        for next_url in _BASE_URLS[1:]:
            results = _retry_sample_data(results, next_url)
            if not results.errors:
                break
        else:
            _handle_final_errors(results)
    return results + already_downloaded
