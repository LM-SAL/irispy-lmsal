0.2.0 (2022-04-14)
==================

Breaking Changes
----------------

- Due to a bug in how older versions of ``irispy-lmsal`` handled the WCS of IRIS files, the readers have been changed.
  This has meant that some of the older API has changed in ways that are not recoverable.
  If your old code has broken, please do get in touch, I want to ensure the old API still works as much as humanly possible.
  Future plans are to restore the older API but that will take some time.

0.1.1 (2022-02-17)
==================

Features
--------

- Added a ``irispy.utils.wobble_movie`` to create a wobble movie. It does need FFMPEG to be installed.

0.1.0 (2022-01-14)
==================

First formal release of ``irispy-lmsal``.

Please note there are parts of this library that are still under going development and will be updated as time
goes on.
There is also a lot of work to be done on the documentation and some of the funciotns in the ``utils`` module
do not function.
