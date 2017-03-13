SimCADO Workshop
================

| University of Vienna
| `Department of Astrophysics <https://astro.univie.ac.at/en/home/>`__
| Türkenschanzstraße 17
| 1170 Vienna
| Austria

| Start: 13:00, 16th February, 2017
| End: 14:00, 17th February, 2017

Summary
-------

A **very** hands-on and informal workshop for using SimCADO to simulate
science images with MICADO. The primary audience are the MICADO and
MAORY science teams, however anybody who wishes to learn the basics of
SimCADO is welcome to attend. However, each participant **MUST** bring
their own science case (and the relevant data to attack it) to the
workshop (more on this below). Python 3 is the language of SimCADO and
therefore also the language of the workshop.

**Please read the `Prerequisits and
Preparation <#prerequisits-and-preparation>`__ section if you plan on
attending**

Format: Hack-a-thon
-------------------

The primary goal of this workshop is to get people comfortable using
SimCADO on their own so that they can continue developing and
investigating their own science cases with SimCADO after the workshop.
As such, the workshop will be run in the style of a hack-a-thon: small
groups programming together and helping each other to get the most out
of the software - similar to the dot-Astronomy and SPIE Hack-Days.

We will begin with a 1-2 hour introduction to the finer points of
SimCADO and work through a series of examples. Then everyone will split
into groups and begin pushing their own science cases through SimCADO.
Oliver and I will move from group to group helping where we can and
giving ideas on how best to attack the different science cases with
SimCADO.

It must be stressed that this is a **very informal** meeting. There is
little in the way of an official programme and coffee breaks can be
taken whenever.

Programme
---------

16th of February:
~~~~~~~~~~~~~~~~~

-  13:00 Welcome and organisational things
-  13:15 Introduction to the finer points of SimCADO
-  13:45 Worked Example 1: Point sources
-  14:15 Worked Example 2: Extended Sources
-  14:45 Coffee!!
-  15:00 Science case presentations
-  15:15 Hack-a-thon begins
-  xx:xx Coffee on the go
-  18:00 “Official” End of Day

17th of February:
~~~~~~~~~~~~~~~~~

-  09:00 Major problems round table
-  10:00 Hack-a-thon continues
-  xx:xx Coffee on the go
-  13:00 Wrap up and suggestions for SimCADO Beta

Venue
-----

| `Department of Astrophysics <https://astro.univie.ac.at/en/home/>`__
| University of Vienna

| Türkenschanzstraße 17,
| 1180 Wien,
| Austria

Prerequisits and Preparation
----------------------------

To participate in this workshop, you will need the following:

-  a laptop (preferably with 8GB of RAM) with:

   -  Python 3
   -  Numpy, Scipy, Astropy, SimCADO

      -  optional: Matplotlib, Poppy, PySynPhot, photutils

-  a science case idea
-  data for your science case

   -  Spatial information. Any one of the following:

      -  e.g. an image (FITS or otherwise) for extended source
      -  or brightness profile and object dimensions
      -  or list of (x,y) coordinates for point sources

   -  Spectral information. Any one of the following:

      -  e.g. Spectrum/Spectra for various objects
      -  or Photometric information: Magnitudes, colours, surface
         brightness
      -  or Fluxes (preferably in ph/s/m2/nm, i.e. F\_lambda)

Please make sure you have a science case and the data to attack it.
While spectators are welcome to attend, the primary reason for this
“work”-shop is to work. Therefore to get the most out of the time,
please **prepare your data and your laptop** before arriving.

A series of tutorials and worked examples will be posted on the SimCADO
website in the coming weeks. It will be advantageous to look through
these before the workshop begins.

Registration
------------

To register, or just to show interest, please send an email to:
kieran.leschinski@univie.ac.at

LOC
---

-  Kieran Leschinski
-  Oliver Czoske
-  Werner Zeilinger
-  Rainer Köhler
-  Michael Mach


Participants and their topics
-----------------------------

-  Ric Davies  
   Looking at observing the Antennae galaxy redshifted to z=1 so that the PaBeta
   filter can be used to detect Halpha
   
-  Suzanne Ramsay  
   Looking into massive star forming regions with stellar spectra surrounded by
   nebulosity
   
-  Davide Massari   
   Relative astrometry: what is the accuracy achievable for a globular cluster 
   science case without any distortion term, but with the only PSF modelling 
   procedure as source of uncertainty?
   
-  Simona Paiano  
   Study the simulation of a QSO (mag 17 in H band) and its 
   fainter elliptical host galaxy (mag 19 or fainter)
   
-  Maximilian Fabricius  
   10 year motions in the core of Omega Cen
   
-  Gijs Verdoes Kleijn
   The moon for calibration purposes
   
-  Giuliana Fiorentino
   Comparisons to MAORY simulations

-  Michele Perna
   Detections of QSOs
   
-  Natascha Förster Schreiber
   Super star clusters in high redshift galaxies


Results of the Workshop
------------------------
-  Ric Davies
   Made a light curve for simulated images of a pulsar with HTR imaging and a 
   windowed readout.  
   Redshifted the Antennae galaxy to z=1 and observed them
   
-  Suzanne Ramsay  
   Did the ground work for making young clusters and the associated nebulocity
   
-  Davide Massari and Giuliana Fiorentino  
   Confirmed that we are on the same level as the MAORY simulations  

-  Maximilian Fabricius
   Created images of the movement of stars in Omega Cen for a period of 10 years.
   Movement is visible, and we should be able to already do astrometry down to 
   10km/s with a 1 year baseline

-  Simona Paiano and Michele Perna
   Looked into resolving the faint host galaxies of quasars
   
-  Natascha Förster Schreiber  
   Built the code base to simulate super star clusters in high redshift
   galaxies. 


Ideas from the Discussion round at the end
-------------------------------------------   
   
Questions
~~~~~~~~~~
* What filters are included? 
* How to use an image (i.e. source_from_image)
    * Diagrammatically explain the scaling

Blackboxes
~~~~~~~~~~~
* what is "lam"
* more examples on the docs
    * definitely for the "must-have" parameters
   

Things to improve for next time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* stars() takes ra,dec coordinates + FOV centre
* rename TC_HAWKI_H"RG
* Start off with examples
* pixel_scale of source_from_image - deal with the grid
* option in simcado.source.galaxy
    * drop cosmologocial info (Oliver) - only use angles and apparent magnitude
* background values - iterate with Ric
* user must specify filter to observe it
* include function to return spectrum for magnitude - e.g. mag_to_spec()
    * or method for Source() to scale based on magnitude, without passing any spectra
* Zeropoints in the header info
* exact pixel poistions in the Source() (e.g <Source>.xpix, .ypix)
* improve noise models 
    * document how to use FPA_PIXEL_MAP
* make sure all default.config parameters are there
* Future: Keywords are the same as the pipeline keywords
    * Look into what are standards from ESO etc
* Include field dependent MAORY SCAO PSFs
* Future: Flat field (Wolfgang)


