# Download

SimCADO can be downloaded from this link 

[http://homepage.univie.ac.at/kieran.leschinski/SimCADO/SimCADO-0.2dev.zip]()

## Installation
​
To install it, download SimCADO from the link above and use the standard `pip` call to install it

`$ pip install --user SimCADO-0.2dev.zip`

Alternatively give the full URL to pip and let it do the downloading for you

`$ pip install --user http://homepage.univie.ac.at/kieran.leschinski/SimCADO/SimCADO-0.2dev.zip`

**Note** that SimCADO will need to download several hundreds of MBs of instrument data into the install directory. Hence why we use the `--user` flag when installing via `pip`. If you want to keep SimCADO in your normal packages directory, then you will need to give python root access while updating SimCADO's data files.

## Getting up-to-date data for SimCADO

SimCADO is being developed along side MICADO. To keep the files that SimCADO uses as fresh as possible, you should run the following command, at the very least **the first time you start SimCADO** 

	>>> import simcado
	>>> simcado.get_extras()

This downloads the latest versions of all the files SimCADO uses behind the scenes, e.g. PSF files, transmission curves etc. If you haven't used SimCADO for several months, chances are there are new data available. It is therefore advantageous to run this command periodically.

## Adding variation to the Detector noise

By default SimCADO only supplies a single H4RG noise frame. Hence if you plan on generating a series of images, all images will have the same *detector noise* pattern. Don't worry, photon shot noise is still completely random, however if you stack 1000s of simulate images, the singular nature of the 



## Dependencies
 
### Needed:
* numpy >1.10.4
* scipy >0.17
* astropy >1.1.2

### Optional (but recommended)
* matplotlib
* poppy >0.4.0
​
