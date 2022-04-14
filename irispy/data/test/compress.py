"""
Short script I used to create the test files in this folder.
"""


def compress(files):
    from astropy.io import fits
    from scipy.ndimage import zoom

    for file in files:
        hdus = fits.open(file)
        for hdu in hdus:
            hdu.verify("fix")
            if "CADPL_DV" in hdu.header:
                print(hdu.header["CADPL_DV"])
                del hdu.header["CADPL_DV"]
            if "CADEX_DV" in hdu.header:
                print(hdu.header["CADEX_DV"])
                del hdu.header["CADEX_DV"]
            if isinstance(hdu, fits.hdu.table.TableHDU):
                continue
            if "NAXIS3" in hdu.header:
                data = []
                for i in range(hdu.data.shape[0]):
                    data.append(zoom(hdu.data[i], 0.1, order=0))
                hdu.data = data
                hdu.header["NAXIS1"] = hdu.data.shape[2]
                hdu.header["NAXIS2"] = hdu.data.shape[1]
                hdu.header["NAXIS3"] = hdu.data.shape[0]
        hdus.writeto(f"{file}.fits")
