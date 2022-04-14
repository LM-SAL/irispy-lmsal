.. _known_issues:

************
Known Issues
************

This page documents commonly known issues, issues here is defined broadly and refers to oddities or specifics of how ``irispy-lmsal`` or the Python ecosystem works that could catch anyone out.

Per-exposure WCS metadata
=========================

The IRIS FITS files contain the per-exposure WCS metadata (reference coordinate and PCij matrix) in an extension, while the primary header has only values averaged over the observation
For example, OBS 4204700138 SJI file, ``CRVAL`` in the primary header is (1.93840, 1.96290), but the reference coordinate
is actually (-3.56638378,  1.69258388), which is retrieved from the appropriate index in the ``XCENIX``/``YCENIX`` arrays in the extension.

This was (and partially still is) a problem for the ``irispy-lmsal`` package, which assumed that the reference coordinate is the same for all exposures.
This has been taken care of by the SJI reader but not fully by the SP reader.

In future, the goal is to create a gWCS to take account of this.

Spectrogram WCS
===============

For different programs that are ran for the slit spectrograph, a sit and start has the CDELT of 0 arcseconds in the X direction.
This is not allowed by the WCS standard, so we set its value to 1e-10 arcsec in the X direction, which essentially fakes the WCS calculation to get it to work if there is no rotation.
However, the PCij matrix used is derived from the SJI, with square pixels, so the PCij matrix is a pure rotation.
This means that one gets the correct answer only if one does the matrix multiplication in the wrong order: first by PCij and then by CDELTs.

We work around this by modifying the PC_ij matrix to have the correct skew.
Since the X CDELT is 1e-10 arcsec, the inverse is thankfully not infinity.
Using equation 187 in `Calabretta & Greisen 2002
<https://www.aanda.org/articles/aa/abs/2002/45/aah3860/aah3860.html>__`, we correct for this.

Note that since these pixels are extremely rectangular, an aspect ratio of ~3e-10, the cross terms in the
PCij matrix are quite small: -3.4e-12 and -3.8e7.
Hopefully, 64-bit floats have enough precision to enable this to work all of the time.
