# Examples with the Source object

The `Source` class is probably the most important class for testing science cases. Therefore spending time on creating accurate `Source` representations of the object of interest is key to getting good results with SimCADO.

Basically a `Source` object represents photon sources using lists of positions (`.x, .y`), a list of unique spectra (`.spectra`) and a list of references which match each photon source to a spectrum in the list of spectra (`.ref`). All sources (extended and point source) can be decomposed into these lists. The advantage of using this approach is that objects with highly similar spectra can both reference the same position in `.spectra`, thereby reducing the number of spectra that need to be manipulated during a simulation. 

## My first Source object
To begin with it is probably easiest to let SimCADO generate a `Source` object. The convenience function `simcado.source.star()` will generate a `Source` object containing a single star. In this case, we'll choose an G2V star with a K-band magnitude of 20, placed 5 arcsec above the centre of the focal plane:

    >>> star_1 = simcado.source.star(mag=20, filter_name="K", spec_type="G2V", position=(5,0))

`star_1` the coordinates of the star are held in the arrays `.x` and `.y`, and the G2V spectrum in `.spectra´.

    >>> star_1.x[0], star_1.y[0]
    (5, 0)

The spectrum for `star_1` is held in `.spectra` and the central wavelength of each of the spectral bins is in `.lam`. 

    >>> star_1.lam, star_1.spectra[0]
    <insert output here>

**Note** that `.lam` is a (1,n) array where as `.spectra` is a (m,n) array where n is the number of bins in the spectra and m is the number of unique spectra in the `Source` object. 

### Combining `Source` objects
Because a `Source` object is just a collection of arrays, it is easy to add many together with the `+` operator:

    >>> star_2 = simcado.source.star(mag=22, filter_name="K", spec_type="A0V", position=(2,-2))
    >>> two_stars = star_1 + star_2
    >>> two_stars.x, two_stars.y
    ((5, 2), (0, -2))

**Note** this a very trivial example (for which SimCADO has a more elegant function: `simcado.source.stars()`), but it serves to illustrate the main idea. The overloaded `+` operator is very useful for combining objects, e.g. forground stars and background galaxies, in order to get a better representation of the sky.

### Point sources
For generating a field of stars, SimCADO offers a series of convenience functions. Please see the docstring or API documentation for more information on how best to use them.
* `simcado.source.star()`
* `simcado.source.stars()`
* `simcado.source.star_grid()`
* `simcado.source.source_1E4_cluster()`


### Using an image as a template for a `Source` object
If we have an extended source that we wish to simulate, e.g. a galaxy, a nebula, etc. we can use the function `simcado.source.source_from_image()`. The image must be a 2D `numpy.ndarray`, but it can come from anywhere, e.g. a FITS file, generated my another function, or even a MS Paint bitmap image. 

    >>> image_1 = astropy.io.fits.getdata("orion.fits")
    >>> image_1
    array([[0.0, ... 0.0],
           ... ,
           [0.0, ... 0.0]])

Here SimCADO takes a the pixel coordinates of the image and converts them to positions on the focal plane. **Note** the user must specify a plate-scale in arcseconds (`pix_res=`) for the image. Each pixel with a value above a certain threshold (default `flux_threshold=0`) will be used in the `Source` object. The coordinates of these pixels are added to the arrays `.x` and `.y`. 

We also need to provide a spectrum for the image. This spectrum is assumed to be the only spectrum for each pixel in the image. The pixel values are then the intensity assigned to that spectrum at that pixel position. 

SimCADO provides the pickles library for stellar spectra. Unfortunately there aren't any built-in galactic spectra yet - for this the user will need to provide their own spectrum.

    >>> lam, spec_1 = simcado.source.SED("G2V", "K", magnitude=20)
    >>> lam, spec_1
    (array([0.7 ... 2.5]), array([0.0 ... 0.0]))

With `image_1`, `lam` and `spec_1` we can now build a `Source` object for an orion-like nebula that has the spectrum of a sun-like star. 

    >>> simcado.source.source_from_image(image_1, lam, spec_1, pix_res=0.004, flux_threshold=0)

While this example is physically unrealistic, it serves the purpose of showing how to build a `Source` object from an image. The user is


### Images with multipe spectra
In reality assigning a single spectrum to an extended object is of limited use. For a `Source` to be realistic is should contain multiple spectra for objects in different locations. The best way to simulate this with SimCADO is to create a `Source` object for each unique group of objects (e.g. old stellar population, star forming regions, AGN, etc) and then combine them into a single `Source` object with the `+` operator.

As a worked example, lets create a "first-order" approximation to a star forming galaxy. The two major components of this source are 
1. the aged stellar population and,
2. the star forming regions. 

In our (very) crude model the aged stellar population can be approxiated by an ellipse with Gaussian light distribution. As M stars make up the majority of this population, we can assign a M0V spectrum to this population.

    >>> from astropy.convolution import Gaussian2DKernel
	>>> from simcado.source import SED
    >>> 
    >>> old_pop = Gaussian2DKernel(128).array[::3,:]
	>>> m0v_spec = SED()

To illustrate (very crudely) the star forming regions we can create a random distribution of elliptical Gaussians using the `astropy` function `Gasussian2DKernel`:








## Making a `Source` object directly


## SimCADO convenience functions
* `simcado.source.stars()`
* `simcado.source.source_1E4_cluster()`
* `simcado.source.SED()`
* `simcado.source.source_from_image()`
* `simcado.source.Source()`