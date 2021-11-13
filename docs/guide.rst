******************************************************
A quick guide to work with IRIS Level 2 data in Python
******************************************************

This guide is intended to help the solar community to start working with IRIS Level 2 data using Python.
It is especially oriented to those who have a limited knowledge of Python and want to start using this language to analyze IRIS data.

This guide will cover:

- reading IRIS Level 2 files
- load the data of any of the spectral windows
- visualize, inspect, and interact with them

In the future more features will arrive:

- a GUI to interact with the data
- analysis tools
- a way to save and load regions of interest
- reading of L3 IRIS data
- Dopplergrams
- Density Diagnostics
- PSF and filter functions


Currently, ``irispy-lmsal`` does not provide a way to download the data from the `IRIS archive <https://iris.lmsal.com/data.html>`__.
We recommend you browse the catalogue using your browser.

This guide uses some "sample data" from the IRIS archive that can be accessed:

.. code-block:: python

  >>> import irispy.data.sample as sample_data  # doctest: +REMOTE_DATA

Once the IRIS Level 2 data has downloaded, the next step is to read, extract, and inspect them.

For now, there are two separate readers available:

.. code-block:: python

  >>> from irispy.io import read_sji_lvl2  # doctest: +REMOTE_DATA
  >>> from irispy.io import read_spectrograph_lvl2  # doctest: +REMOTE_DATA

In future these will be merged into one reader, but for now you need to use the correct one.

The sample data is from this `observation
<https://www.lmsal.com/hek/hcr?cmd=view-event&event-id=ivo%3A%2F%2Fsot.lmsal.com%2FVOEvent%23VOEvent_IRIS_20211001_060925_3683602040_2021-10-01T06%3A09%3A252021-10-01T06%3A09%3A25.xml>`__.

It is a small observation of a small activate region (NOAA 12880) that contains one sunspot.
Let us recover the header of the raster file and show the description of the observation:

.. code-block:: python

    >>> raster = read_spectrograph_lvl2(sample_data.RASTER, uncertainty=False)  # doctest: +REMOTE_DATA

.. note::
    This, by default will load the data into memory.
    You can pass ``memmap=True`` to avoid this, the data array will be a `numpy.memmap` instead.
    Thus, the data are not loaded in the memory of the system, but written in a temporary file.

.. code-block:: python

    >>> raster  # doctest: +REMOTE_DATA
    <ndcube.ndcollection.NDCollection object at ...>
    NDCollection
    ------------
    Cube keys: ('C II 1336', 'Si IV 1394', 'Mg II k 2796')
    Number of Cubes: 3
    Aligned dimensions: [<Quantity 5. pix> <Quantity 16. pix> <Quantity 548. pix>]
    ...

We get some basic information about the raster file from this, what spectral windows were observed
The size of the cube, the wavelength keys as well.

.. note::

  Notice the raster type here, it is a `ndcube.ndcollection.NDCollection` instance.
  In order to allow easier use o multidimensional data sets (which are common in astronomy), the SunPy-affiliated package `ndcube <https://docs.sunpy.org/projects/ndcube/en/stable/>`_  was developed.
  This package combines N-dimensional data with their corresponding word coordinate system (WCS) information
  allowing the data be sliced, visualized and to undergo WCS transforms to alternative coordinate systems.
  ``irispy-lmsal`` is built upon this.

If we want to check the header of the raster, we can do the following:

Let us check the header of this collection, this is stored as a ``meta`` attribute:

.. code-block:: python

    >>> raster["C II 1336"][0].meta  # doctest: +REMOTE_DATA
    <sunraster.instr.iris.IRISSGMeta object at ...>
    IRISMeta
    --------
    Observatory:                IRIS
    Instrument:         SPEC
    Detector:           FUV1
    Spectral Window:    C II 1336
    Spectral Range:             [1331.70275015 1358.28579039] Angstrom
    Date:                       2021-10-01T06:09:25.090
    OBS ID:                     3683602040
    OBS Description:    Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2
    <BLANKLINE>

Note this is not on the main object but each individual element, in this case the spectral window.
While the SJI files contain just one spectral window per file, the raster files have several spectral windows
per file.

As the SJI IRIS Level 2 data are simpler than the raster files, since they have only one spectral window per
file, we will start with this data instead of the raster file above.

We use the following command to read and load the data from a SJI IRIS Level 2 file:

.. code-block:: python

    >>> iris_sji = read_sji_lvl2(sample_data.SJI_1330, uncertainty=False)  # doctest: +REMOTE_DATA
    >>> iris_sji  # doctest: +REMOTE_DATA
    <irispy.sji.IRISMapCube object at ...>
    <BLANKLINE>
    IRISMapCube
    -----------
    Observatory:                 IRIS
    Instrument:                  SJI
    Bandpass:                    1330.0
    Obs. Start:                  2021-10-01T06:09:24.920
    Obs. End:                    2021-10-01T06:11:44.461
    Instance Start:              2021-10-01T06:09:25.020
    Instance End:                2021-10-01T06:11:37.580
    Roll:                        0.000464606 deg
    Total Frames in Obs.:        20
    IRIS Obs. id:                3683602040
    IRIS Obs. Description:       Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2
    Axis Types:                  [('time', 'time', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM', 'custom:CUSTOM'), ('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat'), ('custom:pos.helioprojective.lon', 'custom:pos.helioprojective.lat')]
    Cube dimensions:             [ 20. 548. 555.] pix
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
    <Time object: scale='utc' format='isot' value=2021-10-01T06:09:24.920>

It possible it might be in a ``"DATE_OBS"`` instead.

The exposure times:

.. code-block:: python

    >>> iris_sji.exposure_time   # doctest: +REMOTE_DATA
    <Quantity [0.50031197, 0.50025398, 0.50023699, 0.50024003, 0.50023901,
               0.50028503, 0.50024903, 0.500269  , 0.50026202, 0.500247  ,
               0.50029403, 0.50021601, 0.50028402, 0.50023901, 0.50024903,
               0.50025803, 0.500283  , 0.50029802, 0.50029498, 0.50027299] s>

In most cases, the exposure times are fixed for all scans in a raster.
However, when automatic exposure compensation (AEC) is switched on and there is a very energetic event (e.g. a
flare), IRIS will automatically use a lower exposure time to prevent saturation in the detectors.

If the exposure time varies, you can get the time-dependent exposure times in seconds from the auxiliary
metadata, second to last HDU in the file with the keys ``"EXPTIMEF"`` and ``"EXPTIMEN"``.

To get arrays of timestamps, or exposure times or "xcenix", that information will be in the ``extra_coords``
attribute.

.. code-block:: python

    >>> iris_sji.extra_coords  # doctest: +REMOTE_DATA
    <ndcube.extra_coords.extra_coords.ExtraCoords object at ...>
    ExtraCoords(time (0) None: TimeTableCoordinate ['time'] [None]:
    ['2021-10-01T06:09:25.020' '2021-10-01T06:09:31.990'
     '2021-10-01T06:09:38.950' '2021-10-01T06:09:45.920'
     '2021-10-01T06:09:52.920' '2021-10-01T06:09:59.890'
     '2021-10-01T06:10:06.860' '2021-10-01T06:10:13.860'
     '2021-10-01T06:10:20.830' '2021-10-01T06:10:27.800'
     '2021-10-01T06:10:34.800' '2021-10-01T06:10:41.770'
     '2021-10-01T06:10:48.740' '2021-10-01T06:10:55.740'
     '2021-10-01T06:11:02.700' '2021-10-01T06:11:09.670'
     '2021-10-01T06:11:16.670' '2021-10-01T06:11:23.640'
     '2021-10-01T06:11:30.610' '2021-10-01T06:11:37.580'],
    ...
    exposure time (0) None: QuantityTableCoordinate ['exposure time'] [None]:
    <Quantity [0.50031197, 0.50025398, 0.50023699, 0.50024003, 0.50023901,
               0.50028503, 0.50024903, 0.500269  , 0.50026202, 0.500247  ,
               0.50029403, 0.50021601, 0.50028402, 0.50023901, 0.50024903,
               0.50025803, 0.500283  , 0.50029802, 0.50029498, 0.50027299] s>,
                slit x position (0) None: QuantityTableCoordinate ['slit x position'] [None]:
    <Quantity [258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541,
               258.75      , 270.74543085, 282.74086427, 294.73629541] arcsec pix>,
                slit y position (0) None: QuantityTableCoordinate ['slit y position'] [None]:
    <Quantity [254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75,
               254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75, 254.75,
               254.75, 254.75, 254.75, 254.75] arcsec pix>)

Understanding a Level 2 FITS file
=================================

The structure of the IRIS Level 2 FITS data file is as follows:

The IRIS Level 2 FITS are multi-extension FITS files.
An extension" refers to a part of the file containing self-consistent information.
This information may be, in the general case, a header or its corresponding data.
The first extension is called ``primary`` and its ``extension number`` is 0.

The extensions in an IRIS Level 2 SJI FITS file has the following numbers:

   - ``0``: header and data corresponding to the spectral images observed by the SJI.
   - ``1``: header and auxiliary 31 values from each exposure taken by the SJI in the spectral band of the file.
     It is an array of float values with dimensions :math:`no. images \times 31`.
   - ``2``: header and extra data from each exposure taken by the SJI in the spectral band of the file.
     It is a record array containing 5 string fields for each exposure.
     The values of each field can be access as the key in a dictionary or as an attribute.
     See example in the last code block of this section.

An IRIS Level 2 raster FITS file has the following extensions:

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

The function `astropy.fits.io` shows the information of the extensions contained in the IRIS Level 2 file.
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
    Observation description:  Very large sparse 16-step raster 15x175 16s   Deep x 0.5 Spatial x 2
    <BLANKLINE>
    Extension No. 1 stores data and header of SJI_1330: 1310.00 - 1350.00 AA

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

As the number of spectral windows in a raster file may vary from an observation to another, a good option to access to the last 2 extensions of the IRIS Level 2 file, is to use a negative index:

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

`More information on the Level 2 data can be found in ITN 26. <https://iris.lmsal.com/itn26/iris_level2.html>`__
