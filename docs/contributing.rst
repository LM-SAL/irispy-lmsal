**************************
Contributing to ``irispy``
**************************

``irispy`` is open-source and community-developed, and we are always glad to welcome new contributors and users.
You can contribute in several ways: by providing feedback, reporting bugs, contributing code, and reviewing pull requests.
There is a role for almost any level of engagement.

Providing Feedback
==================

We could always use more voices and opinions in the discussions about ``irispy`` and its development from both users and developers.
There are several ways to make your voice heard.
You can open issues on our `issue tracker`_, comment in the `SunPy chat room`_, or email us directly at ``nabil.freij@gmail.com``.
Whether it be (non)constructive criticism, inquiries about current or future capabilities, or flattering praise, we would love to hear from you.

Reporting Bugs
==============

If you run into unexpected behavior or a bug please report it.
All bugs are raised and stored on our `issue tracker`_.
If you are not sure whether your problem is a bug, a deficiency in functionality, or something else, you can email ``nabil.freij@gmail.com``.
Ideally, we would like a short code example so we can run into the bug on our own machines.

Contributing Code
=================

If you would like to contribute code, it is strongly recommended that you first discuss your aims with the ``irispy`` community.
We strive to be an open and welcoming community for developers of all experience levels.
Discussing your ideas before you start can give you new insights that will make your development easier, lead to a better end product, and reduce the chances of your work being regretfully rejected because of an issue you weren't aware of (e.g., the functionality already exists elsewhere).

In the rest of this section we will go through the steps needed to set up your system so you can contribute code to ``irispy``.
This is done using `git`_ version control software and `GitHub`_, a website that allows you to upload, update, and share code repositories (repos).
If you are new to code development or git and GitHub you can learn more from the following guides:

* `SunPy Newcomers Guide`_
* `GitHub guide`_
* `git guide`_

The principles in the SunPy guides for contributing code and utilizing GitHub and git are exactly the same for ``irispy``, except that we contribute to the irispy repository rather than the ``sunpy`` one.
If you are a more seasoned developer and would like to get further information, you can check out the `sunpy Developers Guide`_.

Before you can contribute code to irispy, you first need to install the development version of ``irispy``.
To find out how, see :ref:`irispy-dev-install`.

Use US English in documentation, comments, and internal names.
Cite scientific references using the :doc:`bibliography` roles so they link to the bibliography.

Releasing
=========

A release is a tag on ``main``; CI builds and publishes it.

#. Preview the changelog with ``towncrier build --draft --version X.Y.Z``, then run ``towncrier build --version X.Y.Z``, which writes the fragments into ``CHANGELOG.rst`` and deletes them, and merge that change through a pull request.
#. Tag the merge commit ``vX.Y.Z`` (``vX.Y.ZbN`` or ``vX.Y.ZrcN`` for a pre-release) and push the tag.
   Never tag a post-release: the CI workflow ignores ``*dev*``, ``*pre*`` and ``*post*`` tags, so ``v0.9.1post1`` made a GitHub release and a Zenodo version but no PyPI release, and ``setuptools_scm`` then reported ``0.9.1.post2.devN`` on ``main``.
   A fix after a release, a documentation fix included, is a new patch release.
#. The tag runs the CI workflow: the tests, the documentation build and ``build_dists``, which builds the sdist and the wheel and tests the installed wheel; the ``publish`` job then uploads both to PyPI through the ``pypi`` environment with trusted publishing.
   Watch the run and check the release on `PyPI <https://pypi.org/project/irispy-lmsal/>`__.
#. Publish a GitHub release from the tag, with the changelog section as its notes.
   The Zenodo integration archives it under the concept DOI `10.5281/zenodo.10443678 <https://doi.org/10.5281/zenodo.10443678>`__; check that the new version appears with the right version number and authors.
#. Read the Docs builds the tag on its own; check that the new version is active, that ``stable`` points at it and that the changelog page shows the release.
#. The conda-forge bot opens a pull request on the ``irispy-lmsal`` feedstock once the release is on PyPI.
   Check the recipe's run requirements and ``python`` floor against ``dependencies`` and ``requires-python`` in ``pyproject.toml`` before merging it.

.. _issue tracker: https://github.com/LM-SAL/irispy/issues
.. _SunPy Newcomers Guide: http://docs.sunpy.org/en/latest/dev_guide/newcomers.html
.. _GitHub: https://github.com/
.. _git: https://git-scm.com/
.. _GitHub guide: https://github.com/git-guides
.. _git guide: https://git-scm.com/book/en/v2/Getting-Started-Git-Basics
.. _sunpy Developers Guide: http://docs.sunpy.org/en/latest/dev_guide
.. _SunPy chat room: https://app.element.io/#/room/#sunpy:openastronomy.org
