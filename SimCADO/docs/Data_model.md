# Data model

**Recent module name changes**
+ Detector.py --> detector.py
+ LightObject --> source.py
+ PSFCube.py --> psf.py
+ OpticalTrain.py --> optics.py
+ SpectralCurve --> spectral.py
+ PlaneEffects --> spatial.py
+ Usercommands.py --> commands.py
+ DefaultsGenerator.py --> defaults.py
+ InputGenerator.py --> input.py
+ SimulationRun.py --> simulation.py



This document should provide a brief overview of the classes available in each module and the interfaces between them.

Here we describe the following classes:
+ `detector.Detector`
+ `detector.Chip`
+ `optics.OpticalTrain`
+ `source.Source`
+ `psf.PSF`
+ `psf.PSFCube`
+ `psf.UserPSFCube`
+ `spectral.TransmissionCurve`
+ `spectral.EmissionCurve`

We also list the functions available in the following modules
+ `spatial`
+ `utils`



## Detector class
The `Detector` class is essentially a container class for the `Chips`. It contains a list of attributes regarding the physical dimensions of the Detector as well as a list of `Chip` instances.

### Attributes
+ cmds : UserCommands
+ layout : astropy.table.Table
+ chips : list
+ oversample : int
+ fpa_res : float
+ exptime : float
+ ndit : int
+ tro : float
+ signal : numpy.ndarray

### Methods

```
__init__(cmds):
	returns <Detector>
```


```
read_out(output=False, filename=None, chips=None, **kwargs):
	if output=True
		returns <astropy.io.fits.HDUList>
	else
		returns None
```

The `Open`method should be a independent function of the Detector class
```
open(filename, **kwargs)
	returns None
```

```
write(filename=None, **kwargs)
	returns None
```


## Chip class
The `Chip` class contains information on each individual chip. It doesn't care where about what is around it. The primary task of the `Chip` class is to create FITS files using the different readout schemes, and to make sure that the signal is sufficiently noisy. "Pure" signal is placed on the chip when `Source.apply_optical_train()` is called. `Chip` then adds the photon shot noise and readout noise before saving it to disk


### Attributes
+ x_cen, y_cen : float
+ x_len, y_len : int
+ x_min, y_min : float
+ x_max, y_max : float
+ naxis1, naxis2 : int
+ pix_res : float
+ id : int
+ array




### Methods
```
__init__(x_cen, y_cen, x_len, y_len, pix_res, id=None)
	returns <Chip>
```

```
add_signal(signal)
	returns None
```

```
add_uniform_background(emission, lam_min, lam_max, output=False)
	returns None
```

```
apply_pixel_map(pixel_map_path=None, dead_pix=None, max_well_depth=1E5)
	returns None
```

```
reset_chip()
	returns None
```

```
read_out(cmds)
	returns numpy.ndarray
```

### Private methods

```
_apply_saturation(arr)
	returns numpy.ndarray
```

```
_read_out_uptheramp(image, dit, tro=1.3, max_byte=2**30)
	returns numpy.ndarray
```

```
read_out_fast(image, dit)
	returns numpy.ndarray
```

```
read_out_superfast(image, dit)
	returns numpy.ndarray
```

Put in a catch for the case that the self.array.shape is larger than the noise frames read in from the noise FITS file. 
E.g if the noise frames are for a 2k chp and we simulate a 4k chip

```
_read_noise_frame(scmds)
	returns numpy.ndarray
```



## OpticalTrain class
The `OpticalTrain` object contains everything needed to fully represent the optical path between the source object and the readout electronics. For those effects which require input data (e.g. `TransmissionCurve`) or an array to function as a convolution kernel (e.g. `PSFCube`), this data is stored in the `OpticalTrain` object. For effects which are best described by an action (e.g. field rotation), `OpticalTrain` contains wrapper functions for the corresponding functions in the `PlaneEffects` module


### Attributes
+ lam : 1D array
+ lam_res : float
+ psf_size : int
+ tc_mirror : SpectralCurve.TransmissionCurve
+ tc_atmo : SpectralCurve.TransmissionCurve
+ tc_source : SpectralCurve.TransmissionCurve
+ ec_mirror : SpectralCurve.EmissionCurve
+ ec_atmo : SpectralCurve.EmissionCurve
+ ph_mirror : SpectralCurve.EmissionCurve
+ ph_atmo : SpectralCurve.EmissionCurve
+ n_ph_mirror : float
+ n_ph_atmo : float
+ cmds : UserCommands.UserCommands
+ detector : Detector.Detector
+ chips : list
+ adc_shifts : list

### Methods
```
__init__(cmds)
	returns <OpticalTrain>
```

The `read` method should be renamed `read_optical_train` and should be turned into an unbound function in the OpticalTrain module

```
read(filename)
	** Not implemented **
	returns <OpticalTrain>
```

```
save(filename)
	** Not implemented **
	returns None
```

```
apply_derotator(arr):
	returns numpy.ndarray
```

```
apply_wind_jitter(arr)
	returns numpy.ndarray
```


### Private methods
```
_make(cmds=None)
	returns None
```

```
_gen_master_tc(tc_keywords=None, preset=None)
	returns SpectralCurve.TransmissionCurve
```

```
_gen_master_psf(psf_type="Airy")
	return PSFCube.AiryPSFCube
```

```
_gen_adc_shifts()
	return list
```

```
_gen_field_rotation_angle()
	** Not implemented **
	returns float
```

```
_gen_telescope_shake()

	returns PSFCUbe.GaussianPSF
```


## Source class
The `Source` class holds two grounds of numpy.ndarrays which describe the spatial (`.x, .y, .ref, .weight`) spectral (`.lam, .spectra`) characteristics of the object in questions.


### Attributes
+ x, y : numpy.ndarray, numpy.ndarray
+ x_orig, y_orig : numpy.ndarray, numpy.ndarray
+ ref : numpy.ndarray
+ weight : numpy.ndarray
+ lam : numpy.ndarray
+ spectra : numpy.ndarray
- pix_res : float
- exptime : float
- area : float
- params : dict
- info : list


### Methods
```
__init__(filename=None, lam=None, spectra=None, x=None, y=None, ref=None, weight=None, **kwargs)
	returns <Source>
```

```
apply_optical_train(opt_train, chips, **kwargs)
	returns None	
```

```
project_onto_chip(image, chip)
	returns None
```

```
image_in_range(psf, lam_min, lam_max, chip=None, **kwargs)
	returns numpy.ndarray
```

```
photons_in_range(lam_min, lam_max, min_bins=10, mask=None)
	returns numpy.ndarray
```

```
photons_in_range(lam_min, lam_max, min_bins=10, mask=None)
	returns numpy.ndarray
```

```
apply_transmission_curve(transmission_curve)
	returns None
```

This method should be turned into an unbound function
```
read(filename)
	returns None
```

```
write(filename)
	returns None
```


### Private Methods

```
_convert_to_photons()
	returns None
```

```
_from_cube(filename, **kwargs)
	returns None
```

```
_from_arrays(lam, spectra, x, y, ref, weight=None)
	returns None
```

```
_from_arrays(lam, spectra, x, y, ref, weight=None)
	returns None
```



## PSF class
The `PSF` class holds the PSF kernel for a certain wavelength

### Attributes

+ size : int
+ shape : (int, int)
+ pix_res : float
+ array : numpy.ndarray
+ info : list


### Methods

```
__init__(size, pix_res)
	returns <PSF>
```

```
set_array(array, threshold=1e-15)
	returns None
```

```
resize(new_size)
	returns None
```

```
resample(new_pix_res)
	returns None
```

```
convolve(kernel)
	returns None
```


## PSFCube class
The `PSFCube ` class holds a list of `PSF` objects

### Attributes

+ lam_bin_centers : list or numpy.ndarray
+ psf_slices : list
+ info : list


### Methods

```
__init__(lam_bin_centers)
	returns <PSFCube>
```

```
resize(new_size)
	returns None
```

```
resample(new_pix_res)
	returns None
```

```
convolve(kernel_list)
	returns None
```

```
export_to_fits(self, filename, clobber=True, **header_info)
	returns None
```

```
nearest(lam)
	returns numpy.ndarray
```



## UserPSFCube class - (PSFCube)
The `UserPSFCube` class will probably be the most used `PSFCube` class because it enables the user to read in their own PSFs from FITS files. As most users will want to use the most accurate PSF, this will mean they will be using the FITS files generated by the MAORY and SCAO teams

### Attributes
+ All the attributes of PSFCube
+ header : astropy.io.fits.Header

### Methods
```
__init__(self, filename, lam_bin_centers)
	returns <UserPSFCube>
```


## TransmissionCurve class
The primary task of `TransmissionCurve` is two hold two arrays - `.lam` and `.val` - which hold the wavelengths and values for a transmission curve.

### Attributes
+ lam : numpy.ndarray
+ val : numpy.ndarray
+ lam_orig : numpy.ndarray
+ lam_orig : numpy.ndarray
+ params : dict
+ info : list


### Methods
```
__init__(**kwargs)
	returns <TransmissionCurve>
```

```
resample(bins, action="average", use_edges=False, min_step=1E-5)
	returns None
```

```
normalize(val=1., mode='integral')
	returns None
```


### Private methods

```
_get_data()
	returns (numpy.ndarray, numpy.ndarray)
```

## EmissionCurve class - (TransmissionCurve)

### Attributes
+ Attributes from TransmissionCurve
+ factor : float

### Methods
```
__init__(**kwargs)
	returns <EmissionCurve>
```

```
resample(bins, action="sum", use_edges=False)
	returns None
```

```
convert_to_photons()
	returns None
```

```
photons_in_range(lam_min, lam_max)
	returns numpy.ndarray
```



## PlaneEffects module
`PlaneEffects` holds a series of functions which can be applies to 2D images. There is no explicit PlaneEffects class.

### Functions
```
tracking(arr, cmds)
	returns numpy.ndarray
```

```
derotator(arr, cmds)
	returns numpy.ndarray	
```

```
wind_jitter(arr, cmds)
	returns numpy.ndarray		
```

```
adc_shift(cmds)
	returns (list, list)
```

### Private Functions
```
_gaussian_dist(x, mu, sig)
	returns numpy.ndarray
```

```
_linear_dist(x)
	returns numpy.ndarray
```

```
_line_blur(arr, shift, kernel="gaussian", angle=0)
	returns numpy.ndarray
```

```
_rotate_blur(arr, angle, kernel="gaussian")
	returns numpy.ndarray
```



## utils module
The utils module contains a list of functions which are (sometimes) used by more than one module. This means that utils should not call other modules, otherwise we'll start an infinite import loop

### Functions
```
msg(cmds, message, level=3)
	returns None
```

```
read_config(config_file)
	returns dict
```

```
update_config(config_file, config_dict)
	returns dict
```

```
unify(x, unit, length=1)
	returns nummpy.ndarray
```

```
parallactic_angle(ha, de, lat=-24.589167)
	returns float
```

```
parallactic_angle_2(ha, de, lat=-24.589167)
	returns float
```

```
moffat(r, alpha, beta)
	returns float
```

```
poissonify(arr)
	returns numpy.ndarray
```

```
atmospheric_refraction(lam, z0=60, temp=0, rel_hum=60, pres=750, lat=-24.5, h=3064)
	returns float or numpy.ndarray
```

```
nearest(arr, val)
	returns int
```

```
poissonify(arr)
	returns numpy.ndarray
```

```
add_keyword(filename, keyword, value, comment="", ext=0)
	returns None
```











