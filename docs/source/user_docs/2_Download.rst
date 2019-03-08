Download
========

The latest stable version of SimCADO can be downloaded from this link

`<http://www.univie.ac.at/simcado/SimCADO.tar.gz>`_

The latest development version of SimCADO can be found on GitHub:

`<https://github.com/astronomyk/SimCADO>`_


Python 3 vs Python 2
--------------------

.. note:: SimCADO has been programmed in Python 3.

While most of the basic functionality for SimCADO will work with Python
2.7, we haven’t tested it properly. For example, the function
``simcado.install_noise_cube()`` only works in Python 3. This isn’t
critical - it just means that if you want to have read noise variations
in your images, you need to use Python 3. However it won’t crash SimCADO
for single images.

.. See the `F <Features>`__ for a list of the "known" issues when
   running in Python 2.7

.. note::
    A side note: Astropy will stop supporting Python 2.7 in 2019 and the
    official End-of-Life for Python 2.7 is 2020, i.e. no more maintainance.
    We are running under the assumption that SimCADO will (hopefully) still be
    around after 2020, hence why we have concentrated our efforts on developing
    in Python 3.

Installation
------------

To install it, download SimCADO from the link above and use the standard
``pip`` call to install it: ::

    $ pip install SimCADO


.. Alternatively give the full URL to pip and let it do the downloading for you: ::

..    $ pip3 install --user http://www.univie.ac.at/simcado/SimCADO.tar.gz


.. Note::
    SimCADO will need to download several hundreds of MBs of instrument data
    into the install directory. Hence why we use the `--user` flag when
    installing via `pip`. If you want to keep SimCADO in your normal packages
    directory, then you will need to give python root access while updating
    SimCADO's data files.

Dependencies
------------

Required

+-----------+-----------+
| Package   | Version   |
+===========+===========+
| numpy     | >1.10.4   |
+-----------+-----------+
| scipy     | >0.17     |
+-----------+-----------+
| astropy   | >1.1.2    |
+-----------+-----------+
| wget      | >3.0      |
+-----------+-----------+
| requests  | >2.0      |
+-----------+-----------+
| synphot   | >0.1      |
+-----------+-----------+
| pyyaml    |           |
+-----------+-----------+

Optional

+--------------+-----------+
| Package      | Version   |
+==============+===========+
| matplotlib   | >1.5.0    |
+--------------+-----------+
| poppy        | >0.4      |
+--------------+-----------+

All dependencies should be automatically installed when typing the `pip` command above.
In case of any error message you can try to install them manually::

   $ pip install numpy scipy astropy wget matplotlib poppy pyyaml synphot requests

We really recomend to install SimCADO in a python distribution like `anaconda`_ or `miniconda`_
as they allow to have full control over your python environment without interfering with your
system instalation. 

.. _anaconda: https://www.anaconda.com/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html

Getting up-to-date data for SimCADO
-----------------------------------

SimCADO is being developed along side MICADO. To keep the files that
SimCADO uses as fresh as possible, you should run the following command,
at the very least **the first time you start SimCADO**::

    >>> import simcado
    >>> simcado.get_extras()

This downloads the latest versions of all the files SimCADO uses behind
the scenes, e.g. PSF files, transmission curves etc. If you haven’t used
SimCADO for several months, chances are there are new data available. It
is therefore advantageous to run this command periodically.

Adding variation to the Detector noise
--------------------------------------

.. note:: Currently only available for Python 3.

By default SimCADO only supplies a single H4RG noise frame. Hence if you
plan on generating a series of images, all images will have the same
*detector noise* pattern. Don’t worry, photon shot noise is still
completely random, however if you stack 1000s of simulated images with
the same detector noise frame, this pattern will show through. This
problem can be avoided by running the following function the first time
you run simcado:

**We recommend running this command in a separate Python session** as it
can take up to 10 minutes depending on your computer.::

    >>> simcado.install_noise_cube(n=25)

SimCADO contains the code to generate unique detector noise images
(`Rauscher 2015 <http://adsabs.harvard.edu/abs/2015PASP..127.1144R>`__),
however creating a 4k detector noise frame takes about 20 seconds. In
order to avoid wasting time by generating noise for each frame every
time SimCADO is run, ``.install_noise_cube(n)`` generates ``n`` noise
frames and saves them in the SimCADO data directory. In future
simulation one of these detector noise frames are picked at random
whenever a ``<Detector>.read_out()`` is called.

Obviously a trade-off has to be made when running
``.install_noise_cube(n)``. The more noise frames available, the less
systematic noise is visible in the read noise of stacked images.
However, the more frames are generated, the longer it takes. A good
solution is to open a separate window and have SimCADO generate frames
in the background.
