.. _install:

************
Installation
************

Release version
===============

``irispy-lmsal`` is part of the wider ecosystem of scientific Python packages for solar physics and therefore a working installation is more about installing the scientific Python ecosystem.

To install the Minifoge Python distribution `download the executable for your Operation System <https://github.com/conda-forge/miniforge#miniforge3>`__.

The reason we choose Miniconda over Anaconda, is mainly due to the size as Anaconda comes with a full install of packages you probably do not need and this way you have more direct control over what has been installed into your Python virtual environment.

Using Miniforge
---------------

To install, launch a system command prompt or the 'Miniforge Prompt' (under Windows).

.. note::

    We strongly recommend using a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

Now to install ``irispy-lmasl`` within the default conda virtual environment:

.. code-block:: console

    $ conda install irispy-lmsal

This will install ``irispy-lmsal``.

Updating
--------
You can update to the latest version by running:

.. code-block:: console

    $ conda update irispy-lmsal

.. _dev_install:

Development version
===================

This section outlines how to install the development version of ``irispy-lmsal``.

Stable Dependencies Install
---------------------------

Create a conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^
The first step is to create a conda environment (let's call it ``irispy-dev``) in which to install the development version of ``irispy-lmsal``.
This will allow you to keep your root environment clean of development packages.
From the command line, type:

.. code-block:: console

    $ conda create -n irispy-dev pip

This creates the ``irispy-dev`` conda environment with just ``pip``.
Next, you must activate that environment, i.e., switch into it.

.. code-block:: console

    $ conda activate irispy-dev

Clone ``irispy-lmsal`` repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we need to clone the `irispy repository`_ from `GitHub`_ into a directory called ``irispy-git``.
From the directory in which you want ``irispy-git`` to reside.
If you want to develop ``irispy-lmsal``, you will want to fork the repository on GitHub and use that URL in the clone step above.

.. code-block:: console

    $ git clone <personal fork URL> irispy-git
    $ cd irispy-git
    $ git remote add upstream git@github.com:LM-SAL/irispy-lmsal.git

Install ``irispy-lmsal``
^^^^^^^^^^^^^^^^^^^^^^^^

Finally, we can install the development version.

.. code-block:: console

    $ cd irispy-git
    $ pip install -e ".[dev]"

You are now e ready to develop ``irispy-lmsal``.

Notice we install no dependencies or use ``conda`` to install this.
The reason for this is that it is simply easier to use ``pip`` to setup development packages.

At times you might need to get the updated changes, to do so:

.. code-block:: console

    $ git remote update -p

From here, you will need to decide if you need to merge changes or rebase changes when you need to contribute the changes back.

.. _irispy repository: https://github.com/LM-SAL/irispy-lmsal
.. _GitHub: https://github.com/
