"""
Short script I used to create the test FITS files in this folder.
"""


def compress(files: list) -> None:
    from astropy.io import fits
    from scipy.ndimage import zoom

    for file in files:
        hdus = fits.open(file)
        for hdu in hdus:
            hdu.verify("fix")
            # TODO; why is this here?
            if "CADPL_DV" in hdu.header:
                print(hdu.header["CADPL_DV"])  # NOQA: T201
                del hdu.header["CADPL_DV"]
            # TODO; why is this here?
            if "CADEX_DV" in hdu.header:
                print(hdu.header["CADEX_DV"])  # NOQA: T201
                del hdu.header["CADEX_DV"]
            if hdu.data is None:
                continue
            # Can't pop out the array, resizing can cause issues
            # So I remove the data and move on.
            if hdu.data.ndim == 1:
                hdu.data = None
                continue
            if hdu.data.ndim == 2:
                factor = (0.1, 1)
            if hdu.data.ndim == 3:
                factor = (0.1, 0.1, 0.1)
            hdu.data = zoom(hdu.data, factor, order=0)
            if hdu.data.ndim == 3:
                hdu.header["NAXIS1"] = hdu.data.shape[2]
                hdu.header["NAXIS2"] = hdu.data.shape[1]
                hdu.header["NAXIS3"] = hdu.data.shape[0]
            if hdu.data.ndim == 2:
                hdu.header["NAXIS1"] = hdu.data.shape[1]
                hdu.header["NAXIS2"] = hdu.data.shape[0]
        hdus.writeto(f"{file.split('.fits')[0]}_test.fits", overwrite=True)
