# Future ideas for SimCADO

This document houses the list of things that need to be implemented or ideas which could be implemented in the future

## To Do List

+ `Detector`
	+ rename the module `detector`
	+ Detector.open move method into an unbound function, returning a `Detector` object
	+ **BUG** - it calls `self.params`which is out-dated
	+ **Potential BUG** - in `Chip._read_noise_frame()` we need to make sure we are prepared for the case that `self.array.shape` > the shape of the noise FITS file

+ `OpticalTrain`	
	+ rename the module `optical_train`
	+ The `read` method should be renamed `read_optical_train` and should be turned into an unbound function in the OpticalTrain module


## Things to be done in the future
+ Unify the names of the read and write methods for each class
+ Change the structure of the `TransmissionCurve` and `EmissionCurve` to replace the `params` dict actual attributes and arguments in `__init__()` E.g. exp_time, pix_res, ...
+ Add a `FilterSet` class to `SpectralClass` that holds all the