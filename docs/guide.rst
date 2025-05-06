.. _guide:

*************************
An ``irispy-lmsal`` guide
*************************

This guide is intended to help the solar community to start working with level 2 data using Python.
It is especially oriented to those who have a limited knowledge of Python and want to start using this language to analyze data.

This guide will cover:

- reading level 2 files
- load the data of any of the spectral windows

In the future more features will arrive:

- a GUI to interact with the data
- a way to save and load regions of interest
- analysis tools
- Density Diagnostics
- Dopplergrams
- PSF and filter functions
- reading of L3 data

Currently, ``irispy-lmsal`` does not provide a way to download the data from the `IRIS archive <https://iris.lmsal.com/data.html>`__.
We recommend you browse the catalogue using your browser.

This guide uses some "sample data" from the IRIS archive that can be accessed:

.. code-block:: python

  >>> import irispy.data.sample as sample_data  # doctest: +REMOTE_DATA

Once the level 2 data has downloaded, the next step is to read, extract, and inspect them.

.. code-block:: python

  >>> from irispy.io import read_files  # doctest: +REMOTE_DATA

The sample data is from this `observation
<https://www.lmsal.com/hek/hcr?cmd=view-event&event-id=ivo%3A%2F%2Fsot.lmsal.com%2FVOEvent%23VOEvent_IRIS_20211001_060925_3683602040_2021-10-01T06%3A09%3A252021-10-01T06%3A09%3A25.xml>`__.

It is a small observation of a small activate region (NOAA 12880) that contains one sunspot.
Let us recover the header of the raster file and show the description of the observation:

.. code-block:: python

    >>> raster = read_files(sample_data.RASTER)  # doctest: +REMOTE_DATA

.. note::
    This, by default will load the data into memory.
    You can pass ``memmap=True`` to avoid this, the data array will be a `numpy.memmap` instead.
    Thus, the data are not loaded in the memory of the system, but written in a temporary file.

.. code-block:: python

    >>> raster  # doctest: +REMOTE_DATA
    <irispy.spectrograph.Collection object at ...>
    <BLANKLINE>
    Collection
    ----------
    Cube keys: (np.str_('C II 1336'), np.str_('Si IV 1394'), np.str_('Mg II k 2796'))
    Number of Cubes: 3
    Aligned dimensions: [5 16 548]
    Aligned physical types: [('meta.obs.sequence',), ...]
    <BLANKLINE>

We get some basic information about the raster file from this, what spectral windows were observed
The size of the cube, the wavelength keys as well.

If we want to check the header of the raster, we can do the following:
Let us check the header of this collection, this is stored as a ``meta`` attribute:

.. code-block:: python

    >>> raster["C II 1336"][0].meta  # doctest: +REMOTE_DATA
    {'SIMPLE': True, 'BITPIX': 16, 'NAXIS': 0, 'EXTEND': True, 'DATE': '2021-11-15', 'COMMENT': "and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H", 'TELESCOP': 'IRIS', 'INSTRUME': 'SPEC', 'DATA_LEV': 2.0, 'LVL_NUM': 2.0, 'VER_RF2': 'L12-2019-08-08', 'DATE_RF2': '2021-11-15T17:16:30.934', 'DATA_SRC': 1.51, 'ORIGIN': 'SDO', 'BLD_VERS': 'V9R41X', 'LUTID': 4.0, 'OBSID': '3683602040', 'OBS_DESC': 'Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2', 'OBSLABEL': '', 'OBSTITLE': '', 'DATE_OBS': '2021-10-01T06:09:25.090', 'DATE_END': '2021-10-01T06:09:51.080', 'STARTOBS': '2021-10-01T06:09:24.920', 'ENDOBS': '2021-10-01T06:11:44.461', 'OBSREP': 5, 'CAMERA': 1, 'STATUS': 'Quicklook', 'BTYPE': 'Intensity', 'BUNIT': 'Corrected DN', 'BSCALE': 0.25, 'BZERO': 7992, 'HLZ': 0, 'SAA': '           0', 'SAT_ROT': 0.00136173, 'AECNOBS': 0, 'AECNRAS': 0, 'DSUN_OBS': 149776000000.0, 'IAECEVFL': 'NO', 'IAECFLAG': 'YES', 'IAECFLFL': 'YES', 'TR_MODE': '', 'FOVY': 182.32, 'FOVX': 14.9658021927, 'XCEN': -321.061, 'YCEN': 390.437, 'SUMSPTRL': 0, 'SUMSPTRN': 2, 'SUMSPTRF': 4, 'SUMSPAT': 2, 'EXPTIME': 0.499401, 'EXPMIN': 0.499328, 'EXPMAX': 0.499448, 'DATAMEAN': 11.3413, 'DATARMS': 21.5989, 'DATAMEDN': 1.51188, 'DATAMIN': -1443.98, 'DATAMAX': 22816.9, 'DATAVALS': 16859200, 'MISSVALS': 1053824, 'NSATPIX': 0, 'NSPIKES': 0, 'TOTVALS': 17913024, 'PERCENTD': 94.117, 'DATASKEW': 22.806, 'DATAKURT': 3317.27, 'DATAP01': -5.68974, 'DATAP10': -2.7028, 'DATAP25': -0.910321, 'DATAP75': 8.3314, 'DATAP90': 28.4133, 'DATAP95': 43.6624, 'DATAP98': 58.6684, 'DATAP99': 67.9847, 'NEXP_PRP': 1, 'NEXP': 16, 'NEXPOBS': 240, 'NRASTERP': 16, 'RASTYPDX': 1, 'RASTYPNX': 1, 'RASRPT': 1, 'RASNRPT': 5, 'STEPS_AV': 0.997720146179, 'STEPS_DV': 5.43998884055e-07, 'STEPT_AV': 1.73733, 'STEPT_DV': 0.100247, 'CADPL_AV': 27.916, 'CADPL_DV': 0.0, 'CADEX_AV': 27.7355, 'CADEX_DV': 0.0116846, 'MISSOBS': 0, 'MISSRAS': 0, 'IPRPVER': 2.54999995232, 'IPRPPDBV': 11.0, 'IPRPDVER': 20130925, 'IPRPBVER': 20211010, 'NWIN': 3, 'TDET1': 'FUV1', 'TDESC1': 'C II 1336', 'TWAVE1': 1335.70996094, 'TWMIN1': 1331.70275015, 'TWMAX1': 1358.28579039, 'TDMEAN1': 0.0460924, 'TDRMS1': 8.75478, 'TDMEDN1': -0.0736134, 'TDMIN1': -697.904, 'TDMAX1': 9416.98, 'TDVALS1': 4271372, 'TMISSV1': 226612, 'TSATPX1': 0, 'TSPIKE1': 0, 'TTOTV1': 4497984, 'TPCTD1': 94.9619, 'TDSKEW1': 86.503, 'TDKURT1': 13743.7, 'TDP01_1': -5.95482, 'TDP10_1': -3.1487, 'TDP25_1': -1.65311, 'TDP75_1': 1.6285, 'TDP90_1': 3.24573, 'TDP95_1': 4.18155, 'TDP98_1': 5.38064, 'TDP99_1': 6.1231, 'TSR1': 1, 'TER1': 548, 'TSC1': 1, 'TEC1': 513, 'IPRPFV1': 100101, 'IPRPGV1': 17, 'IPRPPV1': 100, 'TDET2': 'FUV2', 'TDESC2': 'Si IV 1394', 'TWAVE2': 1393.7800293, 'TWMIN2': 1380.73390787, 'TWMAX2': 1406.73358787, 'TDMEAN2': -0.0916971, 'TDRMS2': 23.7848, 'TDMEDN2': -0.229734, 'TDMIN2': -1443.98, 'TDMAX2': 22816.9, 'TDVALS2': 4204714, 'TMISSV2': 284502, 'TSATPX2': 0, 'TSPIKE2': 0, 'TTOTV2': 4489216, 'TPCTD2': 93.6626, 'TDSKEW2': 196.112, 'TDKURT2': 46462.5, 'TDP01_2': -6.25377, 'TDP10_2': -3.48514, 'TDP25_2': -1.95108, 'TDP75_2': 1.51748, 'TDP90_2': 3.11859, 'TDP95_2': 4.11872, 'TDP98_2': 5.33327, 'TDP99_2': 6.18116, 'TSR2': 1, 'TER2': 548, 'TSC2': 520, 'TEC2': 1031, 'IPRPFV2': 100101, 'IPRPGV2': 17, 'IPRPPV2': 100, 'TDET3': 'NUV', 'TDESC3': 'Mg II k 2796', 'TWAVE3': 2796.19995117, 'TWMIN3': 2783.27137697, 'TWMAX3': 2835.05701538, 'TDMEAN3': 22.8308, 'TDRMS3': 18.7836, 'TDMEDN3': 17.2087, 'TDMIN3': -499.3, 'TDMAX3': 1553.05, 'TDVALS3': 8383114, 'TMISSV3': 542710, 'TSATPX3': 0, 'TSPIKE3': 0, 'TTOTV3': 8925824, 'TPCTD3': 93.9198, 'TDSKEW3': 1.37615, 'TDKURT3': 6.35362, 'TDP01_3': 1.12439, 'TDP10_3': 4.68571, 'TDP25_3': 8.29817, 'TDP75_3': 32.8018, 'TDP90_3': 50.786, 'TDP95_3': 61.1156, 'TDP98_3': 72.91, 'TDP99_3': 80.5779, 'TSR3': 1, 'TER3': 548, 'TSC3': 12, 'TEC3': 1029, 'IPRPFV3': 102689, 'IPRPGV3': 16, 'IPRPPV3': 102, 'KEYWDDOC': 'http://www.lmsal.com/iris_science/irisfitskeywords.pdf', 'HISTORY': 'level2  Version L12-2019-08-08', 'exposure time': <Quantity [0.49932799, 0.499349  , 0.49936101, 0.49937099, 0.49939099,
               0.49937901, 0.49935201, 0.49937299, 0.49937499, 0.49938399,
               0.49937901, 0.499347  , 0.49936   , 0.49937499, 0.49938101,
               0.499387  ] s>, 'exposure FOV center': <SkyCoord (Helioprojective: obstime=None, rsun=695700.0 km, observer=None): (Tx, Ty) in arcsec
        [(-328.04573318, 390.48694024), (-327.04802249, 390.48031793),
         (-326.05034235, 390.47369561), (-325.05269266, 390.46673761),
         (-324.054982  , 390.45901666), (-323.05733235, 390.45138727),
         (-322.05946596, 390.44381891), (-321.06144412, 390.4365252 ),
         (-320.06336037, 390.42917047), (-319.06527704, 390.4217547 ),
         (-318.06722388, 390.41439996), (-317.06920193, 390.40707574),
         (-316.07173934, 390.39917169), (-315.04959247, 390.37637506),
         (-314.05244018, 390.36822686), (-313.05525602, 390.36020074)]>, 'observer radial velocity': <Quantity [-253.12770081, -250.48970032, -247.70889282, -245.26559448,
               -242.44000244, -239.80610657, -237.02940369, -234.58859253,
               -231.76489258, -229.13369751, -226.36210632, -223.87640381,
               -221.10499573, -218.47779846, -215.65640259, -213.2256012 ] m / s>, 'orbital phase': array([0.77429509, 0.77458888, 0.77489805, 0.77517134, 0.77548558,
           0.77577937, 0.77608842, 0.77636158, 0.77667564, 0.77696925,
           0.77727824, 0.77755648, 0.77786553, 0.77815932, 0.77847362,
           0.77874517], dtype='>f8')}


Note this is not on the main object but each individual element, in this case the spectral window.
While the SJI files contain just one spectral window per file, the raster files have several spectral windows per file.

As the SJI level 2 data are simpler than the raster files, since they have only one spectral window per file, we will start with this data instead of the raster file above.

We use the following command to read and load the data from a SJI level 2 file:

.. code-block:: python

    >>> iris_sji = read_files(sample_data.SJI_1330)  # doctest: +REMOTE_DATA
    >>> iris_sji  # doctest: +REMOTE_DATA
    <irispy.sji.SJICube object at ...>
    <BLANKLINE>
    SJICube
    -------
    Observatory:           IRIS
    Instrument:            SJI
    Bandpass:              1330.0
    Obs. Start:            2021-10-01T06:09:24.920
    Obs. End:              2021-10-01T06:11:44.461
    Instance Start:        2021-10-01T06:09:25.020
    Instance End:          2021-10-01T06:11:37.580
    Total Frames in Obs.:  None
    IRIS Obs. id:          3683602040
    IRIS Obs. Description: Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2
    Axis Types:            [('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat', 'time', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM'), ('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat'), ('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat')]
    Roll:                  0.000464606
    Cube dimensions:       (20, 548, 555)
    <BLANKLINE>

Metadata
========

Here we will highlight some of the more important metadata that is available.

We can use it to find out kind of data this is:

.. code-block:: python

    >>> iris_sji.meta["OBS_DESC"]  # doctest: +REMOTE_DATA
    'Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2'

When the observation started:

.. code-block:: python

    >>> iris_sji.meta['STARTOBS']   # doctest: +REMOTE_DATA
    '2021-10-01T06:09:24.920'

It possible it might be in a ``"DATE_OBS"`` instead.

The exposure times:

.. code-block:: python

    >>> iris_sji.exposure_time   # doctest: +REMOTE_DATA
    <Quantity [0.50031197, 0.50025398, 0.50023699, 0.50024003, 0.50023901,
               0.50028503, 0.50024903, 0.500269  , 0.50026202, 0.500247  ,
               0.50029403, 0.50021601, 0.50028402, 0.50023901, 0.50024903,
               0.50025803, 0.500283  , 0.50029802, 0.50029498, 0.50027299] s>

In most cases, the exposure times are fixed for all scans in a raster.
However, when automatic exposure compensation (AEC) is switched on and there is a very energetic event (e.g. a flare), IRIS will automatically use a lower exposure time to prevent saturation in the detectors.

If the exposure time varies, you can get the time-dependent exposure times in seconds from the auxiliary metadata, second to last HDU in the file with the keys ``"EXPTIMEF"`` and ``"EXPTIMEN"``.

To get arrays of timestamps, or exposure times or "xcenix", that information will be in the ``extra_coords`` attribute.

.. code-block:: python

    >>> iris_sji.extra_coords  # doctest: +REMOTE_DATA
    <ndcube.extra_coords.extra_coords.ExtraCoords object at ...>
    ExtraCoords(exposure time (0) None: QuantityTableCoordinate ['exposure time'] [None]:
    <Quantity [0.50031197, 0.50025398, 0.50023699, 0.50024003, 0.50023901,
               0.50028503, 0.50024903, 0.500269  , 0.50026202, 0.500247  ,
               0.50029403, 0.50021601, 0.50028402, 0.50023901, 0.50024903,
               0.50025803, 0.500283  , 0.50029802, 0.50029498, 0.50027299] s>,
                obs_vrix (0) None: QuantityTableCoordinate ['obs_vrix'] [None]:
    <Quantity [-253.13569641, -242.44810486, -231.77319336, -221.11309814,
               -210.41799927, -199.78419495, -189.16329956, -178.50950623,
               -167.91630554, -157.33630371, -146.72239685, -136.17030334,
               -125.63009644, -115.05719757, -104.5714035 ,  -94.14320374,
                -83.69550323,  -73.3214035 ,  -62.97399902,  -52.65399933] m / s>,
                ophaseix (0) None: QuantityTableCoordinate ['ophaseix'] [None]:
    <Quantity [0.77429509, 0.77548558, 0.77667391, 0.77786386, 0.77905941,
               0.78024989, 0.78144038, 0.78263599, 0.78382647, 0.78501666,
               0.78621155, 0.78740203, 0.78859252, 0.78978807, 0.79097688,
               0.79216683, 0.79336196, 0.79455239, 0.79574287, 0.79693335] arcsec>,
                pztx (0) None: QuantityTableCoordinate ['pztx'] [None]:
    <Quantity [-7.97803831e+00, -3.98715830e+00,  3.72256944e-03,
                3.99460268e+00, -7.97803831e+00, -3.98715830e+00,
                3.72256944e-03,  3.99460268e+00, -7.97803831e+00,
               -3.98715830e+00,  3.72256944e-03,  3.99460268e+00,
               -7.97803831e+00, -3.98715830e+00,  3.72256944e-03,
                3.99460268e+00, -7.97803831e+00, -3.98715830e+00,
                3.72256944e-03,  3.99460268e+00] arcsec>,
                pzty (0) None: QuantityTableCoordinate ['pzty'] [None]:
    <Quantity [0.6446346 , 0.66160059, 0.67856681, 0.69553316, 0.6446346 ,
               0.66160059, 0.67856681, 0.69553316, 0.6446346 , 0.66160059,
               0.67856681, 0.69553316, 0.6446346 , 0.66160059, 0.67856681,
               0.69553316, 0.6446346 , 0.66160059, 0.67856681, 0.69553316] arcsec>,
                slit x position (0) None: QuantityTableCoordinate ['slit x position'] [None]:
    <Quantity [258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541] arcsec>,
                slit y position (0) None: QuantityTableCoordinate ['slit y position'] [None]:
    <Quantity [254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75,
               254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75,
               254.75, 254.75, 254.75, 254.75] arcsec>,
                xcenix (0) None: QuantityTableCoordinate ['xcenix'] [None]:
    <Quantity [-321.64163621, -321.64154081, -321.64054553, -321.63951873,
               -321.5924215 , -321.59850309, -321.60135777, -321.56819773,
               -321.55565282, -321.55661478, -321.51550993, -321.5241685 ,
               -321.4984636 , -321.49132346, -321.47172876, -321.48122647,
               -321.46051587, -321.41851219, -321.42161527, -321.42543197] arcsec>,
                ycenix (0) None: QuantityTableCoordinate ['ycenix'] [None]:
    <Quantity [390.41458808, 390.43178122, 390.44696156, 390.46218927,
               390.40669468, 390.41598631, 390.42799954, 390.43424635,
               390.38567211, 390.39919174, 390.41787952, 390.43324879,
               390.40355692, 390.4319302 , 390.43515948, 390.44981385,
               390.41605352, 390.43774154, 390.45774336, 390.47763699] arcsec>)

Understanding a level 2 FITS file
=================================

The structure of the level 2 FITS data file is as follows:

The level 2 FITS are multi-extension FITS files.
An extension" refers to a part of the file containing self-consistent information.
This information may be, in the general case, a header or its corresponding data.
The first extension is called ``primary`` and its ``extension number`` is 0.

The extensions in an level 2 SJI FITS file has the following numbers:

   - ``0``: header and data corresponding to the spectral images observed by the SJI.
   - ``1``: header and auxiliary 31 values from each exposure taken by the SJI in the spectral band of the file.
     It is an array of float values with dimensions :math:`no. images \times 31`.
   - ``2``: header and extra data from each exposure taken by the SJI in the spectral band of the file.
     It is a record array containing 5 string fields for each exposure.
     The values of each field can be access as the key in a dictionary or as an attribute.
     See example in the last code block of this section.

An level 2 raster FITS file has the following extensions:

   -  ``0``: main header with the main information of the observation.
      This header has information about all the spectral windows contained in the file and other relevant and
      general information.
      This extension DOES NOT have spectral data associated with the file.
   -  ``1`` to ``N``: header and data for the N spectral windows contained in the file.
   -  ``N+1``: header and auxiliary 47 values from each exposure considered in the file.
      It is an array of float values with dimensions :math:`no. acquisitions \times 47`.
   -  ``N+2``: header and extra information data from each exposure considered in the file.
      It is a record array containing 9 string fields for each exposure. The values of
      each field can be access as the key in a dictionary or as an attribute.
      See example in the last code block of this section.

The function `astropy.fits.io` shows the information of the extensions contained in the level 2 file.
For a SJI file:

.. code-block:: python

   >>> from astropy.io import fits   # doctest: +REMOTE_DATA
   >>> fits.info(sample_data.SJI_1330)   # doctest: +REMOTE_DATA
    Filename: ...iris_l2_20211001_060925_3683602040_SJI_1330_t000.fits.gz
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     162   (555, 548, 20)   int16 (rescales to float32)
      1                1 ImageHDU        38   (31, 20)   float64
      2                1 TableHDU        33   20R x 5C   [A10, A10, A4, A66, A63]

and for the raster file:

.. code-block:: python

    >>> fits.info("iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits") # doctest: +SKIP
    Filename: iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     215   ()
      1                1 ImageHDU        33   (513, 548, 16)   int16 (rescales to float32)
      2                1 ImageHDU        33   (512, 548, 16)   int16 (rescales to float32)
      3                1 ImageHDU        33   (1018, 548, 16)   int16 (rescales to float32)
      4                1 ImageHDU        54   (47, 16)   float64
      5                1 TableHDU        53   16R x 7C   [A10, A10, A4, A10, A4, A66, A66]

If you would like a bit more information, we have a similar function within `irispy-lmsal`:

.. code-block:: python

    >>> from irispy.io import fitsinfo  # doctest: +REMOTE_DATA
    >>> fitsinfo(sample_data.SJI_1330)  # doctest: +REMOTE_DATA
    Filename: ...iris_l2_20211001_060925_3683602040_SJI_1330_t000.fits.gz
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     162   (555, 548, 20)   int16 (rescales to float32)
      1                1 ImageHDU        38   (31, 20)   float64
      2                1 TableHDU        33   20R x 5C   [A10, A10, A4, A66, A63]

.. code-block:: python

    >>> fitsinfo("iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits") # doctest: +SKIP
    Filename: iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     215   ()
      1                1 ImageHDU        33   (513, 548, 16)   int16 (rescales to float32)
      2                1 ImageHDU        33   (512, 548, 16)   int16 (rescales to float32)
      3                1 ImageHDU        33   (1018, 548, 16)   int16 (rescales to float32)
      4                1 ImageHDU        54   (47, 16)   float64
      5                1 TableHDU        53   16R x 7C   [A10, A10, A4, A10, A4, A66, A66]
    Observation description:  Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2

    Extension No. 1 stores data and header of C II 1336: 1331.70 - 1358.29 AA (FUV)
    Extension No. 2 stores data and header of Si IV 1394: 1380.73 - 1406.73 AA (FUV)
    Extension No. 3 stores data and header of Mg II k 2796: 2783.27 - 2835.06 AA (NUV)

If we now want to recover the main header of any file:

.. code-block:: python

    # The main header of a SJI file
    >>> fits.getheader(sample_data.SJI_1330)  # doctest: +REMOTE_DATA
    SIMPLE  =                    T / Written by IDL:  Mon Nov 15 09:26:15 2021
    BITPIX  =                   16 / Number of bits per data pixel
    NAXIS   =                    3 / Number of data axes
    NAXIS1  =                  555 /
    NAXIS2  =                  548 /
    NAXIS3  =                   20 /
    EXTEND  =                    T / FITS data may contain extensions
    DATE    = '2021-11-15'         / Creation UTC (CCCC-MM-DD) date of FITS header
    COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy
    COMMENT and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H
    TELESCOP= 'IRIS    '           /
    INSTRUME= 'SJI     '           /
    ...

    # The main header of a raster file
    >>> fits.getheader("iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits") # doctest: +SKIP
    SIMPLE  =                    T / Written by IDL:  Mon Nov 15 09:21:38 2021
    BITPIX  =                   16 / Number of bits per data pixel
    NAXIS   =                    0 / Number of data axes
    EXTEND  =                    T / FITS data may contain extensions
    DATE    = '2021-11-15'         / Creation UTC (CCCC-MM-DD) date of FITS header
    COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy
    COMMENT and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H
    TELESCOP= 'IRIS    '           /
    INSTRUME= 'SPEC    '           /
    ...

    # The individual header corresponding to Si IV 1403 in the raster
    >>> fits.getheader("iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits", 2) # doctest: +SKIP
    XTENSION= 'IMAGE   '           / IMAGE extension
    BITPIX  =                   16 / Number of bits per data pixel
    NAXIS   =                    3 / Number of data axes
    NAXIS1  =                  512 /
    NAXIS2  =                  548 /
    NAXIS3  =                   16 /
    PCOUNT  =                    0 / No Group Parameters
    GCOUNT  =                    1 / One Data Group
    ...

The same can be done with the data using `astropy.io.fits.getdata`.

As the number of spectral windows in a raster file may vary from an observation to another, a good option to access to the last 2 extensions of the level 2 file, is to use a negative index:

.. code-block:: python

    # The header corresponding to the extra information extension
    >>> fits.getheader("iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits", -1) # doctest: +SKIP
    XTENSION= 'TABLE   '           / ASCII table extension
    BITPIX  =                    8 / 8 bit bytes
    NAXIS   =                    2 / 2-dimensional ASCII table
    NAXIS1  =                  296 / Number of positions along axis 1
    NAXIS2  =                   16 / Number of positions along axis 2
    PCOUNT  =                    0 / Size of special data area
    GCOUNT  =                    1 / one data group (required keyword)
    TFIELDS =                    7 / Number of fields in each row
    TBCOL1  =                    1 /
    TFORM1  = 'A10     '           /
    TTYPE1  = 'FRMID   '           /
    ...
    # The data for the extra information extension
    >>> data = fits.getdata("iris_l2_20211001_060925_3683602040_raster_t000_r00000.fits", -1) # doctest: +SKIP
    # The names of the records
    >>> data.dtype.names # doctest: +SKIP
    ('FRMID',
     'FUVFDBID',
     'FUVCRSID',
     'NUVFDBID',
     'NUVCRSID',
     'FUVfilename',
     'NUVfilename',
     'FUVtemp',
     'NUVtemp')

We can access to the values of the variables stored in the data corresponding to the extra information extension as an attribute or as a key:

.. code-block:: python

    # An example is the record: "FUVfilename"
    >>> data_extra.FUVfilename # doctest: +SKIP
    chararray(['/irisa/data/level1/2021/10/01/H0600/iris20211001_06092534_fuv.fits',
              '/irisa/data/level1/2021/10/01/H0600/iris20211001_06092706_fuv.fits',
              ...
              '/irisa/data/level1/2021/10/01/H0600/iris20211001_06094981_fuv.fits',
              '/irisa/data/level1/2021/10/01/H0600/iris20211001_06095140_fuv.fits'],
              dtype='<U66')

`More information on the level 2 data can be found in ITN 26. <https://iris.lmsal.com/itn26/iris_level2.html>`__
