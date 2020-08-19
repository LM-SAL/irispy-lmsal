************
Installation
************

.. _dev_install:

Installing the development version
==================================

This section outlines how to install the development version of irispy.
The two primary packages on which irispy relies are `ndcube`_ and `sunraster`_.
Both of these have stable released versions that work with irispy.
However, some developers may want to use the latest updates of these packages in their work on irispy.
Below we will first outline how to install irispy with its stable dependencies, and then with the development versions of ndcube and sunraster.

To install these packages we will use a combination of conda, conda environments, pip and git.
We will assume these are all installed on your current system.
If you do not have anaconda installed, see the `anaconda website`_ for instructions.
The main goal of using anaconda is due to its easy virtual environment support but you can use any Python virtual environment.

If you want to develop irispy, we also suggest forking the repository on GitLab and using that in the clone steps below.

Stable Dependencies Install
---------------------------

Create a conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^
The first step is to create a conda environment (let's call it ``irispy-dev``) in which to install the development version of irispy.
This will allow you to keep your root environment clean of development packages.
From the command line, type:

.. code-block:: console

    conda config --append channels conda-forge
    conda create -n irispy-dev pip

The first line opens a conda channel so that irispy and its dependencies can be installed.
The second line creates the ``irispy-dev`` conda environment with a list of dependencies.
Next, you must activate that environment, i.e. switch into it.

.. code-block:: console

    conda activate irispy-dev

Clone irispy Repository
^^^^^^^^^^^^^^^^^^^^^^^
The second step is to clone the `irispy repository`_ from `GitLab`_ into a directory called ``irispy-git``.
From the directory in which you want ``irispy-git`` to reside, type:

.. code-block:: console

    git clone https://gitlab.com/LMSAL_HUB/iris_hub/irispy-lmsal irispy-git

Install irispy
^^^^^^^^^^^^^^
Finally, we can install the irispy development version.

.. code-block:: console

    cd irispy-git
    pip install -e .

You should now be ready to use irispy.
To check it's installed, open an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

    import irispy

To make sure you have the latest updates, regularly do

.. code-block:: console

    git pull origin master

Development Dependencies Install
--------------------------------

We installed the stable versions of sunraster and ndcube above in order to get get all their dependencies.
Now that is done, the second step is to remove the stable versions of sunraster and ndcube and install the development versions.
CAUTION: Make sure you are in (have activated) the ``irispy-dev`` conda environment otherwise the next steps will remove sunraster and ndcube from the wrong environment.

Clone development versions of sunraster and ndcube
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let's make a directory and then clone (download) the development versions of `sunraster`_ and `ndcube`_ into subdirectories.
Let's call them ``sunraster-git`` and ``ndcube-git``.
As above, if you want to create bugfixes you will also want to fork ndcube and sunraster.

On the command line from the directory in which you want your repos to live, type:

.. code-block:: console

    git clone https://github.com/sunpy/ndcube.git ndcube-git
    git clone https://github.com/sunpy/sunraster.git sunraster-git

If you already have these repos cloned, make sure they are up-to-date but by pulling the latest version of the master branches.
For example, for sunraster, do:

.. code-block:: console

    cd ~/sunraster-git
    git pull origin master

assuming that ``origin`` is the remote pointing to the main sunraster repo, i.e. https://github.com/sunpy/sunpy.git.
The same should be done for ndcube.

Install the Development Versions of SunPy, ndcube and irispy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally you will want to change into those locations and run

.. code-block:: console

    pip install -e .

in each folder to install the current development version of these packages.

.. _SunPy: http://sunpy.org
.. _anaconda website: https://docs.anaconda.com/anaconda/install.html
.. _irispy repository: https://gitlab.com/LMSAL_HUB/iris_hub/irispy-lmsal/
.. _GitLab: https://gitlab.com/
.. _sunraster: https://github.com/sunpy/sunraster
.. _ndcube: https://github.com/sunpy/ndcube
