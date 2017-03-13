Testing SimCADO
================

As per a talk on unit testing that I watched on YouTube, the testing regime 
should tell a story. In essence, if all the code was deleted, we should be able
to reconstruct the functionality of SimCADO based solely on the tests. 

Here I aim to outline how that would work by breaking down the SimCADO framework
into its component units, and then assigning tests to them. The tests will 
include the inteded functionality as well as the boundary cases. Where possible
true negative tests will also be included, i.e. tests which check the proper
response of a unit when the wrong input is given.


The SimCADO story
------------------

Running a simulation with SimCADO involves executing the following steps. We
will use these as the base for structuring the tests. 

#. Create a Source object
#. Create a dictionary of UserCommands
#. Create an OpticalTrain based on the UserCommands
#. Create an Detector array based on the UserCommands
#. Apply the effects of the OpticalTrain on the Source
#. Project the resulting image of the Source into the Detector array
#. Read the Detector images out into a FITS file

Now lets break each of these test categories down:


Create a Source object
-----------------------

The Source object
~~~~~~~~~~~~~~~~~~
The Source object should describe the source of light on the sky. The 
functionality of this object should include:

#. Describe the spatial distribution of the sources on sky in terms of angular 
   separation wrt to the centre of the FoV
#. Describe the spectral energy distribution of the sources on sky in ph/s/m2
#. Connect the spatial information with the spectral information
#. Be able to shift the positions on sky
#. Be able to rotate the positions on sky around an arbitary point
#. Be able to scale the spectral information
#. Be able to scale the object spatialy on sky
#. Be able to redshift the object in the spectral domain

The Source object should have the following properties

#. Object should be combinably with the "+" operator
#. A TransmissionCurve should be applyable to the spectum with the "*" operator
#. Object should be saveable to the hard disk
#. Saved Source objects should be readable


Functions with return Source objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#. A function to return a Source object for a single star
#. A function to return a Source object for a group of stars
#. A function to return a Source object representing a cluster of stars with 
   a mass distribution following that of a regular IMF
#. A function to return a Source object for a galaxy with a spatial distribution
   following a Sersic profile of index negative
#. A function to return a Source object for an empty region of the sky
#. A function to return an array describing the wavelength bins and the spectral
   energy distribution of a common celestial object (star, galaxy)
#. A spectrum to return a Source object described by an intensity profile given 
   in an image object (FITS image, jpeg image, numpy array) and a spectrum
#. A function to return a flat spectrum scaled to a certain magnitude



Create a dictionary of UserCommands
------------------------------------
The UserCommands object should describe everything that happens during the 
simulation. This includes describing how the OpticalTrain and the Detector
should be contructed and how the "observation" should take place. 

The UserCommands should contain the following functionality:

#. Should contain a dictionary with all the Keyword-Value pairs described in 
   a default config file
#. Should be able to read files with sub-sets of the the full command dictionary
#. Should allow updates to the command dictionary
#. Allow the user to save the current command dictionary to file
#. Be able to read a saved command dictionary from file
#. Should check that all file names are existing files in any of the accepted 
   path directories

The module should also contain the functionality to:

#. Write configuration files to disk containing all the necessary information
   to run an observation simulation. Such files will contain information 
   regarding the:
   
   #. control commands
   #. layout of the chips in the detector focal plane array
   #. mirror properties in the telescopes optical path


Create an OpticalTrain based on the UserCommands
-------------------------------------------------

The OpticalTrain class is responsible for describing the path that light takes 
through between the emitive source and the detector chip in the focal plane. 
This means that the OpticalTrain must accurately describe the Atmosphere, the
series of mirrors in the telescope's optical system, the effect of the optical
surfaces within the instrument and the characteristics of the detectors in the 
focal plane. 

The following functionality must be included in any object which attempts to 
describe this interplay between photons and the optical path:

Description of the major parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#. Contain a list of mirrors, their dimensions and temperatures for the telescope
#. Contain a list of mirrors, their dimensions and temperatures for the AO system
#. Contain the number of mirrors in the instrument
#. Contain a list of other optical surfaces in the instrument
#. Contain transmission curves/reflectivity curves for each of the optical 
   surfaces
#. Contain a post-AO PSF, ideally one that varies with wavelength
#. Contain a model of the atmopsheric dispersion
#. 


Purely spectral functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. contain a spectral energy distribution for each major source of background
   light in ph/s
#. be able to alter these SEDs to account for the transmission loses along the
   path to the Detector
#. 



Create an Detector array based on the UserCommands
---------------------------------------------------


Apply the effects of the OpticalTrain on the Source
----------------------------------------------------


Project the resulting image of the Source into the Detector array
------------------------------------------------------------------


Read the Detector images out into a FITS file
----------------------------------------------
































