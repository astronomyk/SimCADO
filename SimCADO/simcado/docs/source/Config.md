# SimCADO configuration for MICADO

SimCADO was made for optical train defined by the E-ELT and MICADO. As such the base configuration is for this instrument/telescope combination.

Although SimCADO allows the user to configure the optical train as they like, it is *highly* recommended to leave the default parameters alone if the goal of the simulations is comparability between scientific use cases.

This document lists the current configuration which best describes SimCADO. Each section describes a physical effect included by SimCADO and **gives the Keyword/Value pair** contained in the configuration file **for MICADO**.

## MICADO optical elements
### Point spread functions


### Atmospheric Dispersion Correction 

`INST_ADC_PERFORMANCE = 100`

`INST_ADC_NUM_SURFACES = 8`

`INST_ADC_TC = "default"`

The number of surfaces in the ADC is set to 8. A transmission coefficient of 0.95 is assumed for each of these surfaces. 

SimCADO implements an inverse ADC. It assumes that the ADC corrects perfectly for the dispersion induced by the atmosphere. If the ADC isn't working as it should, i.e. `INST_ADC_PERFORMANCE` less than 100%, SimCADO 


### Derotator

`INST_DEROT_PROFILE = linear`

`INST_DEROT_PERFORMANCE = 100`


### MICADO mirrors

`INST_NUM_MIRRORS = 8`
	
The coating for the MICADO mirrors is assumed to be similar to the coating used for the E-ELT's mirrors. Thus the same reflectivity curve is used for 

#### Number of reflective surfaces
List of reflective surfaces for wide field imaging mode (4mas/pixel)
* 
* 
* 
List of reflective surfaces for wide field imaging mode (4mas/pixel)
*
*
*



### Cryostat window




### Instrumental distorion

`INST_DISTORTION_MAP = none`


## MICADO Detector effects

### MICADO Detector read modes




## E-ELT Spatial Effects

### Wind jitter and telescope vibrations



## E-ELT Transmission Effects

### E-ELT optical surfaces
