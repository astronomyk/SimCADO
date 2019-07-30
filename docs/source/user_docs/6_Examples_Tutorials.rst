Examples and Tutorials
======================

Here are a list of scripts andnotebooks detailing how to use SimCADO. If you would
like to add a notebook to this collection, please email it to `Kieran
Leschinski`_.

.. _Kieran Leschinski: kieran.leschinski@univie.ac.at

Cheat Sheet
-----------

.. hint:: 

    If you don't like sifting through documentation, try looking through the 
    cheat sheet at some common commands and examples. If you find what you need, 
    then you'll know exactly what to look for in the
    :doc:`../reference/simcado`

    * `PDF version`_ or 
    * `Jupyter Notebook`_


.. _PDF version: ../_static/pdfs/SimCADO_cheatsheet.pdf
.. _Jupyter Notebook: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/SimCADO-cheat-sheet.ipynb



Tutorials
---------

-  `my\_first\_sim.ipynb`_

   An introductory notebook to SimCADO. Topics include: first steps with
   SimCADO, creating :class:`.Source` objects and customising simulations.

.. _my\_first\_sim.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/my_first_sim.ipynb

- `stars_anisocado.ipynb`_ 

    A example of how stars looks like with a field changing PSF as expected in the SCAO mode. Here,
    a bunch of stars are observed at 0, 7, 14 and 21 arcsec off-axis. Also available as a stand-alone script (`stars_anisocado.py`_)

.. _stars_anisocado.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/stars_anisocado.ipynb

.. _stars_anisocado.py: https://github.com/astronomyk/SimCADO/blob/master/docs/source/_static/scripts/stars_anisocado.py

- `star_cluster_fullarray.ipynb`_

    This example simulates a star cluster filling the whole MICADO array with a MAORY-like PSF in the J-band. For that a new J-band PSF was generated with a Strehl of 7%. This PSF can be downloaded from the `data directory`_ (as well as others in H and Ks). Also available as stand-alone script 
    (`star_cluster_fullarray.py`_)

.. _star_cluster_fullarray.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/star_cluster_fullarray.ipynb

.. _star_cluster_fullarray.py: https://github.com/astronomyk/SimCADO/blob/master/docs/source/_static/scripts/star_cluster_fullarray.py

.. _data directory: https://github.com/astronomyk/SimCADO/tree/master/data

.. star_cluster_fullarray.py: https://github.com/astronomyk/SimCADO/blob/master/docs/source/_static/scripts/star_cluster_fullarray.py

- `extended_source.ipynb`_ 

    We put here a local galaxy and see how it would look like with MICADO at z=2. We use a SCAO PSF. 
    Stand-alone script: `extended_source.py`_

.. _extended_source.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/extended_source.ipynb

.. _extended_source.py: https://github.com/astronomyk/SimCADO/blob/master/docs/source/_static/scripts/extended_source.py

- `galaxy_globclusters.ipynb`_

    We simulate here a disk galaxy surrounded by (proto-) globular clusters and using a number of observing possibilites with SimCADO. Stand-alone script:  `galaxy_globclusters.py`_

.. _galaxy_globclusters.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/galaxy_globclusters.ipynb

.. _galaxy_globclusters.py: https://github.com/astronomyk/SimCADO/blob/master/docs/source/_static/scripts/galaxy_globclusters.py   


   
Science Cases
-------------

-  `Limiting_Magnitudes.ipynb`_

   Examples for how to use ``limiting_mags()`` from the sub-module 
   ``simcado.simulation``. 

.. _Limiting_Magnitudes.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/Limiting_Magnitudes.ipynb

-  `Betelgeuse.ipynb`_

   A quick look at whether it would be possible to observe Betelgeuse (J
   = -3 mag) with SimCADO without destroying the detectors. Short
   answer: yes.
   
.. _Betelgeuse.ipynb:     https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/Betelgeuse.ipynb   
   
-  `StarDestroyer.ipynb`_

   If General Tarkin were to attack Earth and has hidding in orbit
   around the moon in his spaceship, could the E-ELT see it, and would
   we we know what it was?

.. _StarDestroyer.ipynb:  https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/StarDestroyer.ipynb
   
-  `Examples_Galaxies.ipynb`_

   A quick look to some of the cababilities of SimCADO to work with extended sources

.. _Examples_Galaxies.ipynb:  https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/Examples_Galaxies.ipynb

   
PSF experiments
------------------

- `SimCADO\_PSF\_examples.ipynb`_

  A notebook showing the different PSFs available in SimCADO and doing a bit of analysis of their characteristics.



.. _SimCADO\_PSF\_examples.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/SimCADO\_PSF\_examples.ipynb 

-  `E-ELT\_2-Phase\_Mirror.ipynb`_

   As the primary mirror of the E-ELT was innitially thought to be built in 2 phases, this
   notebook uses SimCADO to have a look at what that means for the
   resulting E-ELT PSF.

.. _E-ELT\_2-Phase\_Mirror.ipynb: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/POPPY_EELT.ipynb   
   
.. note::
    the PSFs generated here are only approximations based on a diffraction limited core combined with a Seeing Halo. These are meant only as a guide.

.. -  `Filter\_Wings.ipynb`_
..
..    A experiment to investigate the transmission coefficients that the
..   filters in MICADO will need to have in order to restrict flux leakage
..   through the filter wings to less than 1% of the flux coming through
..   the filter.
.. 
..   .. _Filter\_Wings.ipynb:  https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/Filter_Wings.ipynb   

Workshop notebooks
------------------

1. `Setting up SimCADO`_

   Some basics about running SimCADO
   
.. _Setting up SimCADO:   https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/1_Setting_up_SimCADO.ipynb   
   
2. `Working with point sources`_

   The basics for making and combining point source objects
   
.. _Working with point sources: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/2_Working_with_Point_Sources.ipynb   
   
3. `Notes on simulation with sub-pixel resolution`_

   A couple of examples on how to get SimCADO to do sub-pixel scale simulations

   .. _Notes on simulation with sub-pixel resolution: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/4_Sub-pixel_shifting.ipynb   
   
   
4. `Working with extended sources`_

   Creating a galaxy field with a large galaxy (the Antennae) and adding a series of background galaxies
   Here are the `FITS file used in the simulation`_

.. _Working with extended sources: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/Antennae_galaxy_field.ipynb 
.. _FITS file used in the simulation: ../_static/downloads/Antennae_data.zip

