# Future ideas for SimCADO

This document houses the list of things that need to be implemented or ideas which could be implemented in the future

## To Do List

+ `UserCommand`
    + add the use of a default UserCommands to all classes that need it

+ `Source`
    + add a return array option to `apply_optical_train()` - really to all methods
    + overload the + oberator to combine two Source objects

+ `Detector`
	+ rename the module `detector` - Done
	+ Detector.open move method into an unbound function, returning a `Detector` object
	+ **BUG** - it calls `self.params`which is out-dated
	+ **Potential BUG** - in `Chip._read_noise_frame()` we need to make sure we are prepared for the case that `self.array.shape` > the shape of the noise FITS file
    + Add the ability to have a different sensitivity frame for each chip

+ `OpticalTrain`	
	+ rename the module `optics` - Done
	+ The `read` method should be renamed `read_optical_train` and should be turned into an unbound function in the OpticalTrain module
    + define and create "pre-compiled" optical train files
    + The ADC is not the only thing dictating the spacing between wavelength slices - LOOK INTO THIS!

## Things to be done in the future
+ Unify the names of the read and write methods for each class
+ Change the structure of the `TransmissionCurve` and `EmissionCurve` to replace the `params` dict actual attributes and arguments in `__init__()` E.g. exp_time, pix_res, ...
+ Add a `FilterSet` class to `SpectralClass` that holds all the