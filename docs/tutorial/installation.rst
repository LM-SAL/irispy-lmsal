.. _irispy-tutorial-installing:

************
Installation
************

This is the first chapter in the ``irispy`` tutorial, and by the end of it you should have a working installation of Python and ``irispy``.
For further information and alternative methods for installing ``irispy`` beyond the recommended approach outlined below, refer to sunpy's documentation (:ref:`sunpy-topic-guide-installing`).

.. warning::

    Please be aware that the package name on pypi and conda-forge is `irispy-lmsal <https://pypi.org/project/irispy-lmsal/>`__ to avoid name clashes with other packages.
    However, the package is imported as ``irispy`` and is referred to as ``irispy`` in the documentation.

Installing Python
=================

There are many ways to install Python, but even if you have Python installed somewhere on your computer we recommend following these instructions anyway.
That's because we will create a new Python virtual environment.
As well as containing a Python installation, this virtual environment provides an isolated place to install Python packages (like ``irispy``) without affecting any other current Python installation.
If you already have Python and ``conda`` working you can skip the next section.

`If you are using Anaconda, we recommend that you uninstall it as the default package channel(s) have a restrictive license which means you might not be able to use it for free <https://sunpy.org/posts/2024/2024-08-09-anaconda/>`__.
Instead, we recommend that you use miniforge which is a minimal installer that set ups ``conda`` with the ``conda-forge`` channel, which is free to use for everyone.
If you are using miniforge, you can skip the next section.

.. _irispy-tutorial-installing-miniforge:

Installing miniforge (and conda)
================================

If you don't already have a Python installation then we recommend installing Python with `miniforge <https://github.com/conda-forge/miniforge/#miniforge>`__.
Miniforge will install ``conda`` and automatically configure the default channel (a channel is a remote software repository) to be ``conda-forge``, which is where ``irispy`` is available.

First, download the installer for your system and architecture from the links below:

.. grid:: 3

    .. grid-item-card:: Linux

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh>`__

        `aarch64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh>`__

        `ppc64le <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-ppc64le.sh>`__

    .. grid-item-card:: Windows
        :link: https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>`__

    .. grid-item-card:: Mac

        `arm64 (Apple
        Silicon) <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh>`__

        `x86-64 <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh>`__

Then select your platform to install miniforge:

.. tab-set::

    .. tab-item:: Linux & Mac
        :sync: platform

        For Linux & Mac, run the downloaded script above using the following command:
        ``bash <filename>``. The following should work:

        .. code-block:: console

            bash Miniforge3-$(uname)-$(uname -m).sh

        Once the installer has completed, restart your terminal or log-out and log-in if the changes don't take effect.

    .. tab-item:: Windows
        :sync: platform

        Double click the executable file downloaded from the links above.

        Once the installer has completed you should have a new "Miniforge Prompt" entry in your start menu.

In a new terminal (Miniforge Prompt on Windows) run ``conda list`` to test that the install has worked.

Installing irispy
=================

To install ``irispy``, start by launching a terminal (under a UNIX-like system) or the Miniforge Prompt (under Windows).
Now we will create and activate a new virtual environment to install ``irispy`` into:

.. code-block:: bash

    $ conda create --name irispy
    $ conda activate irispy

In this case the virtual environment is named 'irispy'.
Feel free to change this to a different environment name.

The benefit of using a virtual environment is that it allows you to install packages without affecting any other Python installations or versions on your system.
This also means you can work on multiple projects (research or coding) with different package requirements without them interfering with each other.

.. dropdown:: Click here if you haven't installed miniforge
    :color: warning

    If you have installed miniforge or are using Anaconda you need to configure conda to get your packages from conda-forge as well as the defaults channel.

    You should no longer use the defaults channel at all, see `this blog post <https://sunpy.org/posts/2024/2024-08-09-anaconda/>`__ for details as to why.
    Therefore, if you are using Anaconda or miniconda we would suggest you uninstall it and install miniforge in its place.

    We also appreciate this isn't going to be possible for everyone, so what follows is our best instructions for how to proceed if you are using miniconda or Anaconda.

    The commands you need to run to add conda-forge and make it the default location to install conda packages from are:

    .. code-block:: bash

        $ conda config --add channels conda-forge
        $ conda config --set channel_priority strict

    These commands are taken from the
    `conda-forge documentation <https://conda-forge.org/docs/user/introduction/#how-can-i-install-packages-from-conda-forge>`__.

    Running these commands affect all the environments in your conda installation, critically, including the base Anaconda environment.
    We highly recommend that you do not install new packages, upgrade packages or use your base environment.
    Instead create new environments for all your projects, as you are much less likely to run into any pitfalls while using `multiple channels <https://conda-forge.org/docs/user/tipsandtricks/#multiple-channels>`__ by doing this.

Now that we have a fresh virtual environment, we can proceed with installing ``irispy``:

.. code-block:: bash

    $ conda install irispy-lmsal

This will install ``irispy`` and all of its dependencies.
If you are planning on using irispy in jupyter notebooks we also recommend you install the ``ipywidgets`` and ``itables`` packages.

To ensure that ``irispy`` was installed correctly, run the following command:

.. code-block:: bash

    $ conda list irispy-lmsal

This checks if ``irispy`` was installed correctly.

If you want to install another package later, you can run ``conda install <package_name>``.

.. _irispy-install-pip:

Installing with pip
-------------------

``irispy`` is also on PyPI, so in any Python environment (Python 3.13 or newer) you can run:

.. code-block:: bash

    $ pip install irispy-lmsal

.. _irispy-install-extras:

Optional dependencies
---------------------

Some functions need packages that are not installed by default.
They are grouped in extras, which you add to the package name, for example ``pip install "irispy-lmsal[cosmic-rays,rgb]"``.
Each of these functions raises an error that names the missing package and the extra.
Not every package is available on conda-forge; the ``rsliding`` package, the default backend of ``remove_cosmic_rays``, is only on PyPI, so install it with ``pip`` even in a conda environment.

.. list-table::
   :widths: 20 30 50
   :header-rows: 1

   * - Extra
     - Packages
     - What it unlocks
   * - ``cosmic-rays``
     - ``rsliding``, ``astroscrappy``
     - `~irispy.utils.cosmic_rays.remove_cosmic_rays` with the ``"rsliding"`` (default) or ``"astroscrappy"`` method
   * - ``rgb``
     - ``colorsynth``
     - The false-color images of `irispy.utils.rgb`
   * - ``all``
     - The two extras above
     - Everything above
   * - ``tests``, ``docs``, ``dev``
     - Test, documentation and all of the above
     - Running the tests and building the documentation, see :ref:`irispy-dev-install`

.. _irispy-dev-install:

Installing a development version
================================

To work on the code, fork `the repository <https://github.com/LM-SAL/irispy>`__ on GitHub, clone your fork and install it in editable mode into a new environment (Python 3.13 or newer):

.. code-block:: bash

    $ git clone https://github.com/<your-username>/irispy.git
    $ cd irispy
    $ conda create --name irispy-dev python=3.14 pip
    $ conda activate irispy-dev
    $ python -m pip install -e ".[dev]"
    $ python -m pip install pre-commit tox

The ``dev`` extra installs the test and documentation dependencies and all optional dependencies.
Run the tests with ``pytest``, the style checks with ``pre-commit run --all-files`` and the documentation build with ``tox -e build_docs``, and add a changelog entry in ``changelog/`` for a user-visible change.

Now we've got a working installation of ``irispy``, in the next few chapters we'll look at some of the basic data structures ``irispy`` uses for representing IRIS SJI and SG data.
