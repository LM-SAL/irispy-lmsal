.. _install:

************
Installation
************

Release version
===============
``irispy-lmsal`` is part of the wider ecosystem of scientific Python packages for solar physics and therefore a working installation is more about installing the scientific Python ecosystem.
If you do not currently have a working scientific Python distribution this guide will set you up with the Miniconda, which makes it easy to install and manage your scientific Python packages.

To install the Miniconda Python distribution follow the `instructions <https://docs.conda.io/en/latest/miniconda.html>`__.
Although Miniconda makes it simple to switch between Python versions, we recommend that new users install the latest Python 3.x version of Miniconda.

The reason we choose Miniconda over Anaconda, is mainly due to the size as Anaconda comes with a full install of packages you probably do not need and this way you have more direct control over what has been installed into your Python virtual environment.
Furthermore, you bypass the need for the conda resolver to sort out your root environment which should make conda faster to use.

Using Miniconda
---------------
To install, launch a system command prompt or the 'Anaconda Prompt' (under Windows).
First configure conda for to add the `conda-forge channel <https://conda-forge.org/>`__::

    conda config --add channels conda-forge
    conda config --set channel_priority strict

and now to install ``irispy-lmasl`` within the default conda virtual environment::

    $ conda install irispy-lmsal

This will install ``irispy-lmsal`` and every package it needs to function.

.. note::
    We strongly recommend using a `Python virtual environment <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`__ or a `conda virtual environment. <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__

Updating
--------
You can update to the latest version by running::

    conda update irispy-lmsal

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

    conda create -n irispy-dev pip

This creates the ``irispy-dev`` conda environment with just ``pip``.
Next, you must activate that environment, i.e., switch into it.

.. code-block:: console

    conda activate irispy-dev

Clone ``irispy-lmsal`` repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now we need to clone the `irispy repository`_ from `GitLab`_ into a directory called ``irispy-git``.
From the directory in which you want ``irispy-git`` to reside, type:

.. code-block:: console

    git clone https://gitlab.com/LMSAL_HUB/iris_hub/irispy-lmsal irispy-git

If you want to develop ``irispy-lmsal``, we suggest forking the repository on GitLab and using that in the clone step above.

Install ``irispy-lmsal``
^^^^^^^^^^^^^^^^^^^^^^^^
Finally, we can install the development version.

.. code-block:: console

    cd irispy-git
    pip install -e ".[dev]"

You should now be ready to use ``irispy-lmsal``.
To check it's installed, open an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

    import irispy

To make sure you have the latest updates, regularly do

.. code-block:: console

    git pull origin master

Development dependencies
------------------------
We installed the stable versions of many packages.
If you want to install development versions of any package you can do the following steps:

- Git clone the source code of the package into a directory called ``package-name-git``.
  e.g., ``git clone https://github.com/sunpy/sunraster.git sunraster-git``
- Change into the directory ``package-name-git``.
  e.g., ``cd sunraster-git``
- Install it using ``pip install -e .``.

.. _irispy repository: https://gitlab.com/LMSAL_HUB/iris_hub/irispy-lmsal/
.. _GitLab: https://gitlab.com/
