# Getting Started with SimCADO
SimCADO can be super easy to use, or super complicated. The level of complexity is completely up to the user. Regardless of your intended use for SimCADO, it's probably a good idea to at least have a vague understanding of what is going on during a simulation.

## Behind the scenes of SimCADO
SimCADO uses 4 main classes during a simulation: 

* `Source` holds spatial and spectral information about the astronomical source of photons, e.g. galaxy, star cluster, etc.
* `OpticalTrain` contains information on the various elements along the optical path, e.g. mirrors reflectivity curves, PSFs, instrumental distortion, etc.
* `Detector` represents the focal plane detector array and contains information on the electronic characteristics of the detector chips and their physical positions.
* `UserCommands` is a dictionary of all the important keywords needed by SimCADO to run the simultationm, e.g. `OBS_EXPTIME` (exposure time) or `INST_FILTER_TC` (filter curve)


## Source
For full details, please see the [API](API) and examples of [Source Objects](examples/Source)

The `Source` class is probably the most important class for testing science cases. Therefore spending time on creating accurate `Source` representations of the object of interest is key to getting good results with SimCADO. `Source` objects can be created from scratch, with functions provided by SimCADO, or by loading in a pre-existing `Source`-FITS file.

**Note: SimCADO is CaSe SensITIVe!** SimCADO has the class `simcado.Source()` and the module `simcado.source`. These should not be confused. `simcado.source` is the module which contains the class `Source` and all the helper functions for creating various types of `Source` objects. In fact the source code for the  class `Source` is actually at `simcado.source.Source`, however to make things easy, `Source` is available directly as`simcado.Source()`. Be careful and remember `simcado.Source != simcado.source`.

For a description of the `Source` object, and the `source` module, see [How SimCADO works](deep_stuff/SimCADO/#source).

### Loading a pre-existing `Source` object

To load in a pre-existing `Source`, specify the keyword `filename=` when initialising the `Source` object.

```
>>> import simcado as sim
>>> my_src = sim.Source(filename="star_grid.fits")
```

`Source`-FITS files have a very specific file format, so it's best to only import files that were generated directly from other `Source` objects. It's a chicken/egg scenario, which is why the next section deals with creating `Source` objects in memory. For a description of the file format for saved `Source` objects, see ["File Format of saved Source objects"](deep_stuff/SimCADO/#source).

### Making a `Source` with SimCADO's in-built functions

The `simcado.source` module provides an ever-increasing series of functions to create `Source` objects in memory. These include, (from `simcado.source`)

* `.empty_sky()`
* `.star(mag, filter_name="K", ...)`
* `.stars(mags, x, y, ...)`
* `.star_grid(n, mag_min, mag_max, ...)`
* `.source_1E4_Msun_cluster(distance=50000, ...)`
* `.source_from_image(images, lam, spectra, pix_res, ...)`

Two useful functions here are `.stars()` and `.source_from_image()`

* `stars()` takes a list of magnitudes (and optionally spectral types) and positions for a common broad-band filter (default is "K") and generates a `Source` object with those stars in the field.

```
>>> x, y = [-2.5, 0.7, 16.3], [3.3, -0.2, 25.1]
>>> mags, spec_types = [25,21,28], ["K0V", "A0III", "G2V"]
>>> filt = "H"
>>>
>>> my_src = sim.source.stars(mags=mags, x=x, y=y, filter_name=filt, spec_types=spec_types)
```

* `source_from_image()` creates a `Source` based on a 2D numpy array provided by the user. The 2D array can come from anywhere, e.g. the data from a FITS image, a BITMAP image, from memory, etc. Alongside the image, the user must provide a spectrum (plus a vector with the bin centres) and the pixel field of view (e.g. 0.004 arcsec for MICADO). SimCADO then extracts all pixels from the image which have values above `flux_threshold` (defualt is 0) and saves these pixel coordinates. The spectrum provided is then connected to these pixel, and scaled by the pixel value.

```
>>> # ... Create an image - a circle with a radius of 20 pixels on a grid 200 pixel wide
>>> XX = np.array([np.arange(-100,101)]*201) 
>>> im = np.sqrt(XX**2 + XX.transpose()**2)
>>> im[im>20] = 0; im[im>0] = 1
>>>
>>> # ... Pull in the spectrum for a G2V star with K=20
>>> lam, spec = simcado.source.SED("G2V", filter_name="K", magnitude=20)
>>>
>>> # ... Make the source object
>>> my_src = sim.source.source_from_image(images=im, lam=lam, spectra=spec, pix_res=0.004)
```

`source_from_image()` can also take a list of images if different spectra are to be assigned to each image. An example of this maybe for galaxies. The older population might be represented by image with an ellipse on it, while the positions of the star forming regions are shown on a different image with random scattered blobs. In this case, both images can be passed in a list to `images` and the array passed to `spectra` must have dimensions (2,n) where n is the length of the spectra. **Note** the spectra must be on the same grid and be the same length.


### Creating a `Source` object from scratch
To create a `Source` object from scratch, we initialise the object by passing 5 (or 6) arrays. All the parameter names must be specified.

`
sim.Source(lam=, spectra=, x=, y=, ref=, [weight=])
`

where: 
+ `x, y` - [each a `numpy.ndarray`]. Coordinates for each point source in the image in units of [arcsec] from the focal plane centre

+ `lam` - [`numpy.ndarray`]. An array with the centre of the wavelength bins in [um] for each unique spectrum

+ `spectra` - [`numpy.ndarray`]. An (n, m) array holding n spectra, each with m values. Default units are [ph/s]
Note - `lam` and `spectra` should use a constant bin width. Variable bin widths leads to unpredictable results.

+ `ref` - [`numpy.ndarray`]. An array to connect the point source at `x[i]`, `y[i]` to a unique spectrum at `spectra[j]`, i.e. `ref[i] = j`

Optional keywords can be specified:

+ `weight` - [`numpy.ndarray`], optional. If two sources share the same spectrum, but are at different distances or have different luminosities a scaling factor can be specified to the spectrum when applied to each specific point source.
+ `units` [default `"ph/s"`] is the units for the spectra, i.e. n phontons per second per spectral bin. The size of the spectral bins is resolution of the `.lam` array.


### Combining two (or more) `Source` objects
`Source` objects can be created in different ways, but the underlying table-structure is the same. Therefore adding `Source` objects together means simply combining tables. The mathematical operator `+` can be used to do this:

```
>>> # ... create a A0V star at (0,0) and a G2V star at (5,-5)
>>> star_A0V = sim.source.star(20, position=(0,0), spec_type="A0V")
>>> star_G2V = sim.source.star(20, position=(5,-5), spec_type="G2V")
>>> 
>>> src_combi = star_A0V + star_G2V
>>> 
>>> print(src_combi.x, src_combi.y)
[0 5] [ 0 -5]
```

By adding different `Source` objects together, it is possible to build up complex objects that will be representative of the observed sky, e.g. old + new galaxy stellar population + gas emission + foreground stars

See [examples](examples/Source) for how to use the `*` and `-` operators with a `Source` object


### Saving a `Source` object to disk
The `Source` object is saved as a FITS file with two extensions. See [How SimCADO works](in_depth/SimCADO) for more on the file structure.

```
>>> src_combi.write("my_src.fits")
```
The file can be read in at a later time by specifying `filename=` when initialising a `Source` object - as stated above

```
>>> my_src = sim.Source(filename="my_src.fits")
```


### In-built `Source` object for a 10<sup>4</sup> M<sub>O</sub> cluster
As a test object, SimCADO provides the function, with all distances in parsecs:

`sim.source.source_1E4_Msun_cluster(distance=50000, half_light_radius=1)`


## OpticalTrain
## Detector
## UserCommands

## Simulating with SimCADO
## The quick, the dirty and the ugly
### Creating `Source` objects

## Building your own simulation run
### Creating `UserCommand` objects
### Creating `OpticalTrain` objects
### Creating `Detector` array objects
### Running the simulation

## Things to watch out for





