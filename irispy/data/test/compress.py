"""
Short script I used to create the test FITS files in this folder.

Observation is https://www.lmsal.com/hek/hcr?cmd=view-event&event-id=ivo%3A%2F%2Fsot.lmsal.com%2FVOEvent%23VOEvent_IRIS_20210905_001833_3620258102_2021-09-05T00%3A18%3A332021-09-05T00%3A18%3A33.xml

Original details
----------------

SnS at AR 12863
OBS 3620258102: Medium sit-and-stare

Where
-----
x,y:    -53",-400"
Max FOV:        60"x60"
Target: AR

Raster
------
FOV:    0"x60"
Steps:  1872x0"
Step Cad:       9.3s
Raster Cad:     9s, 1 ras
Linelist:       v36_01

SJI
---
FOV:    60"x65"
1330:   34s, 520 imgs
1400:   28s, 624 imgs
2796:   28s, 624 imgs
2832:   167s, 104 imgs
"""


def compress(files: list) -> None:
    from scipy.ndimage import zoom  # NOQA: PLC0415

    from astropy.io import fits  # NOQA: PLC0415

    for file in files:
        hdus = fits.open(file)
        for hdu in hdus:
            hdu.verify("fix")
            # TODO: why is this here?
            if "CADPL_DV" in hdu.header:
                print(hdu.header["CADPL_DV"])  # NOQA: T201
                del hdu.header["CADPL_DV"]
            # TODO: why is this here?
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
            elif hdu.data.ndim == 3:
                factor = (0.1, 0.1, 0.1)
            hdu.data = zoom(hdu.data, factor, order=0)
            if hdu.data.ndim == 2:
                hdu.header["NAXIS1"] = hdu.data.shape[1]
                hdu.header["NAXIS2"] = hdu.data.shape[0]
            elif hdu.data.ndim == 3:
                hdu.header["NAXIS1"] = hdu.data.shape[2]
                hdu.header["NAXIS2"] = hdu.data.shape[1]
                hdu.header["NAXIS3"] = hdu.data.shape[0]
        hdus.writeto(f"{file.split('.fits')[0]}_test.fits", overwrite=True)
