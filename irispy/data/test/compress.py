"""
Creates the test files in this folder.
"""


def compress_sji(files):
    from astropy.io import fits
    from scipy.ndimage import zoom

    for file in files:
        hdus = fits.open(file)
        for hdu in hdus:
            hdu.verify("fix")
            if isinstance(hdu, fits.hdu.image.PrimaryHDU):
                if "NAXIS3" in hdu.header:
                    hdu.data = zoom(hdu.data, 0.1, order=0)
                    hdu.header["NAXIS1"] = hdu.data.shape[2]
                    hdu.header["NAXIS2"] = hdu.data.shape[1]
                    hdu.header["NAXIS3"] = hdu.data.shape[0]
                    hdu.header["CDELT1"] = str(float(hdu.header["CDELT1"]) * (1 / 0.1))
                    hdu.header["CDELT2"] = str(float(hdu.header["CDELT2"]) * (1 / 0.1))
                    hdu.header["CDELT3"] = str(float(hdu.header["CDELT3"]) * (1 / 0.1))
        hdus.writeto(f"{file}.fits", overwrite=True)


def compress_raster(files):
    from astropy.io import fits
    from scipy.ndimage import zoom

    for file in files:
        hdus = fits.open(file)
        for hdu in hdus:
            hdu.verify("fix")
            if isinstance(hdu, fits.hdu.image.ImageHDU):
                if "NAXIS3" in hdu.header:
                    hdu.data = zoom(hdu.data, 0.1, order=0)
                    hdu.header["NAXIS1"] = hdu.data.shape[2]
                    hdu.header["NAXIS2"] = hdu.data.shape[1]
                    hdu.header["NAXIS3"] = hdu.data.shape[0]
                    hdu.header["CDELT1"] = str(float(hdu.header["CDELT1"]) * (1 / 0.1))
                    hdu.header["CDELT2"] = str(float(hdu.header["CDELT2"]) * (1 / 0.1))
                    hdu.header["CDELT3"] = str(float(hdu.header["CDELT3"]) * (1 / 0.1))
        hdus.writeto(f"{file}.fits", overwrite=True)


if __name__ == "__main__":
    import glob

    files = glob.glob("/home/nabil/Data/IRIS/*SJI*.fits")
    compress_sji(files)
    files = glob.glob("/home/nabil/Data/IRIS/*raster*.fits")
    compress_raster(files)
