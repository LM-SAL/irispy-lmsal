=========
Changelog
=========

.. note::

    This README was adapted from the pytest changelog readme under the terms of the MIT license.

This directory contains "news fragments" which are short files that contain a small **ReST**-formatted text that will be added to the next ``CHANGELOG``.

The ``CHANGELOG`` will be read by users, so this description should be aimed at users instead of describing internal changes which are only relevant to the developers.

Make sure to use full sentences with correct case and punctuation, for example::

    Added ``irispy.utils.moments.average_window``, which averages a spectrogram cube over a wavelength window.

Name APIs as ``double-backtick`` literals, not Sphinx links, so that an entry still builds after the API it names is removed.
Describe only what changed for users since the last release: a fix or change to something not yet released belongs in the fragment that added it, ``breaking`` and ``removal`` are only for released APIs, and tests, CI and other changes users never see need no entry.

Each file should be named like ``<PULL REQUEST>.<TYPE>[.<COUNTER>].rst``, where ``<PULL REQUEST>`` is a pull request number, ``COUNTER`` is an optional number if a PR needs multiple entries with the same type and ``<TYPE>`` is one of:

* ``breaking``: A change which requires users to change code and is not backwards compatible. (Not to be used for removal of deprecated features.)
* ``feature``: New user facing features and any new behavior.
* ``bugfix``: Fixes a reported bug.
* ``doc``: Documentation addition or improvement, like rewording an entire session or adding missing docs.
* ``deprecation``: Feature deprecation
* ``removal``: Feature removal.
* ``trivial``: A packaging change users can see, such as the package metadata or the files in the wheel.

So for example: ``123.feature.rst``, ``456.bugfix.rst``.

If you are unsure what pull request type to use, don't hesitate to ask in your PR.

Note that the ``towncrier`` tool will automatically reflow your text, so it will work best if you stick to a single paragraph, but multiple sentences are OK.
You can install ``towncrier`` and then run ``towncrier build --draft`` if you want to get a preview of how your change will look in the final release notes.
