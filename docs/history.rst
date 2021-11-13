
************************
``irispy-lmsal`` history
************************

``IRISpy`` was originally developed as a community effort by members of the SunPy team.
The ``IRISpy`` package was never released to the public and exists only in a development form.
Recent changes within the way that SunPy is organized have seen the more general functionality of ``IRISpy``
superseded in the SunPy package by a more instrument agnostic module for analysis of imaging spectrometer data
named `sunraster <https://github.com/sunpy/sunraster/>`__.

``IRISpy`` has therefore been taken over by the IRIS team at LMSAL and renamed ``irispy-lmsal``.
As it exists today, ``irispy-lmsal`` has two main functions, which both do essentially the same thing, but for slit jaw imager (SJI) and Spectrograph respectively:

* Those being functions to load in fits files extract the data from them into the format of an `ndcube <https://docs.sunpy.org/projects/ndcube/en/stable/introduction.html>`__ which
  the preferred method of looking at NDdata with Sunpy.

An alternate way of looking at this is that ``irispy-lmsal`` will form the instrument specific portion of ``sunraster`` and thus provide the gateway for pulling IRIS data into that package.
However, the line between packages hasn't been really drawn formally: the feeling right now is that anything that is instrument/satellite/mission independent should be in ``sunraster``, but these things cannot be separated that easily, and these discussions are likely to continue.

``irispy-lmsal`` therefore remains in an unreleased format, and is under heavy development, and as such, great
care should be taken when using these codes.
However, we expect rapid development and the addition of new features on a continuous basis.
