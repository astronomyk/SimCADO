Notebooks
=========

Here are a list of notebooks detailing how to use SimCADO. If you would
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
.. _Jupyter Notebook: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/SimCADO-cheat-sheet.ipynb



Tutorials
---------

-  `my\_first\_sim.ipynb`_

   An introductory notebook to SimCADO. Topics include: first steps with
   SimCADO, creating :class:`.Source` objects and customising simulations.

.. _my\_first\_sim.ipynb: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/my_first_sim.ipynb
   
Science Cases
-------------

-  `Limiting_Magnitudes.ipynb`_

   Examples for how to use ``limiting_mags()`` from the sub-module 
   ``simcado.simulation``. 

.. _Limiting_Magnitudes.ipynb: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/Limiting_Magnitudes.ipynb

-  `Betelgeuse.ipynb`_

   A quick look at whether it would be possible to observe Betelgeuse (J
   = -3 mag) with SimCADO without destroying the detectors. Short
   answer: yes.
   
.. _Betelgeuse.ipynb:     http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/Betelgeuse.ipynb   
   
-  `StarDestroyer.ipynb`_

   If General Tarkin were to attack Earth and has hidding in orbit
   around the moon in his spaceship, could the E-ELT see it, and would
   we we know what it was?

.. _StarDestroyer.ipynb:  http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/StarDestroyer.ipynb
   
   
MICADO experiments
------------------

-  `E-ELT\_2-Phase\_Mirror.ipynb`_

   As the primary mirror of the E-ELT may be built in 2 phases, this
   notebook uses SimCADO to have a look at what that means for the
   resulting E-ELT PSF.

.. _E-ELT\_2-Phase\_Mirror.ipynb: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/POPPY_EELT.ipynb   
   
.. note::
    the PSFs generated here are only approximations based on a diffraction limited core combined with a Seeing Halo. These are meant only as a guide.

-  `Filter\_Wings.ipynb`_

   A experiment to investigate the transmission coefficients that the
   filters in MICADO will need to have in order to restrict flux leakage
   through the filter wings to less than 1% of the flux coming through
   the filter.
   
.. _Filter\_Wings.ipynb:  http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/Filter_Wings.ipynb   

Workshop notebooks
------------------

1. `Setting up SimCADO`_

   Some basics about running SimCADO
   
.. _Setting up SimCADO:   http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/1_Setting_up_SimCADO.ipynb   
   
2. `Working with point sources`_

   The basics for making and combining point source objects
   
.. _Working with point sources: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/2_Working_with_Point_Sources.ipynb   
   
3. `Notes on simulation with sub-pixel resolution`_

   A couple of examples on how to get SimCADO to do sub-pixel scale simulations
   
.. _Notes on simulation with sub-pixel resolution: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/4_Sub-pixel_shifting.ipynb   
   
   
4. `Working with extended sources`_

   Creating a galaxy field with a large galaxy (the Antennae) and adding a series of background galaxies
   Here are the `FITS file used in the simulation`_

.. _Working with extended sources: http://nbviewer.jupyter.org/url/www.univie.ac.at/simcado/_static/downloads/Antennae_galaxy_field.ipynb 
.. _FITS file used in the simulation: ../_static/downloads/Antennae_data.zip

