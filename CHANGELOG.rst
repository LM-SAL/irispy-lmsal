0.1.4 (2022-09-26)
==================

Features
--------

- Added a timestamp to each frame of the wobble movie.
  You will need to set the ``timestamp`` keyword to be `True`.
- Added a ``wobble_cadence`` keyword to override the default wobble cadence of 180 seconds.

0.1.3 (2022-05-22)
==================

Features
--------

- Added V5 and V6 support for ``get_iris_response``. It also does not download the files anymore.

Breaking Changes
----------------

- API of ``get_iris_response`` has changed:
  ``pre_launch`` has gone, use ``response_version=2`` instead.
  ``response_file`` keyword has been removed, it will use files provided by the package instead.
  ``force_download`` was removed as the function now does not download any files.

0.1.2 (2022-05-02)
==================

Features
--------

- Tweaked ``irispy.utils.wobble_movie`` to remove limits on the metadata.
- Pin ``sunraster`` version due to Python version incompatibilities.

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
