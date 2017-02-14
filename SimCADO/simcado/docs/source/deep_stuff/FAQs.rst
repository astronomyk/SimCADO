SimCADO FAQs
============

Here are some answers to known issues with SimCADO.

Work around for failing :func:`~simcado.detector.install_noise_cube` with Python 2.7
-------------------------------------------------------------------------------------

The problem lies with Python 2.7. The noise cube code is 3rd party code
that only works on Python 3 and I haven’t had a chance to dig into that
code yet to find the problem

Option 1
~~~~~~~~

Generate a noise cube in Python 3. First install python3
::

    $ pip3 install <path-to-simcado.zip> 
    $ python3
    
Then create the noise cube

::

    >>> import simcado
    >>> sim.detector.make_noise_cube(num_layers=25, filename='FPA_noise.fits', multicore=True)

Option 2
~~~~~~~~

Download a 15 slice noise cube from my Google Drive folder and save it
in the simcado/data folder.
`https://drive.google.com/file/d/0B8SnQxFuNeltVVc1RTJ5ZFBDQ0k/view?usp=sharing 
<https://drive.google.com/file/d/0B8SnQxFuNeltVVc1RTJ5ZFBDQ0k/view?usp=sharing>`__

.. note::
    Your simcado/data folder can be found by printing the ``__pkg_dir__``
    variable: 
    ::
    
        >>> simcado.utils.__pkg_dir__

Copy the new noise cube into the Python 2.7 simcado/data folder.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default SimCADO looks for the noise cube in its data directory -
``<Your Python 2.7 Directory>/lib/python/site-packages/simcado/data/``

::

    $ cp ./FPA_noise.fits  <Your Python 2.7 Directory>/lib/python/site-packages/simcado/data/FPA_noise.fits

No access to the ``simcado/data`` folder?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you can’t save files into the simcado/data directory (or you can’t be
bothered finding it), you can use the FPA\_NOISE\_PATH keyword when
running a simulation to point simcado to the new noise cube file

::

    >>> simcado.run(my_source, ..... , FPA_NOISE_PATH="<path/to/new>/FPA_noise.fits")

or if you’re using a UserCommands object to control the simulation:

::

    >>> cmds = simcado.UserCommands()
    >>> cmds["FPA_NOISE_PATH"] = "<path/to/new>/FPA_noise.fits"

Surface brightness scaling in SimCADO
-------------------------------------

.. admonition:: Question

    Is it true that if the array that I pass to simcado has a value of 1,
    then simcado will simulate exactly a SB = ``m`` mag/sq.arcsec pixel?
    What happens if the passed numpy array has a different pixel size. Are
    the counts in each pixel then scaled according to their different
    surface area?

    
To answer the question on scaling we need look at docstring for
``source_from_image``:

::

    source_from_image(images, lam, spectra, plate_scale, oversample=1,
                      units="ph/s/m2", flux_threshold=0,
                      center_pixel_offset=(0, 0),
                      conserve_flux=True)

The 3 parameters of interest here are ``plate_scale``, ``oversample``
and ``conserve_flux``.

-  ``plate_scale`` is the plate scale of the image. So if we take a
   HAWK-I image, we have a plate\_scale of 0.106 arcsec. SimCADO needs
   this so that it can work out how many photons are coming in per pixel
   in the image you provide.

-  ``oversample``. If oversample stays at 1 (default), then SimCADO will
   generate a “point source” for each pixel - i.e. every 0.106 arcsec
   there will be one source emitting with the intensity of that pixel.
   This is sub-optimal if we want to use an extended object, as SimCADO
   will turn the image into a grid of point sources with a spacing equal
   to ``plate_scale``. Therefore we need to over sample the image so
   that SimCADO makes at least 1 light source per pixel. Hence to get
   down to 4mas for the SimCADO wide field mode, we would need to
   oversample the HAWK-I image by a factor of 0.106/0.004 = 26.5. I.e.
   ``oversample=26.5``.

-  ``conserve_flux``. If this is ``True`` (default), then when SimCADO
   oversamples the image, it multiplies the new image by a scaling
   factor
   ``= sum(orig_image) / sum(new_image). Thus a pixel with the value 1 in the original image will now have a value``\ =(1
   / 26.5)^2\ ``. If``\ conserve\_flux=False\`, this scaling factor is
   not applied. **Note:** The scaling doesn’t affect the spectrum
   associated with the image at all.

Sub-pixel accuracy with SimCADO
-------------------------------

There are two ways to go about getting SimCADO to simulate at the
sub-pixel level (Nov 2016, only one is working so far):

-  Turn on SimCADO’s oversample mode with the keyword SIM\_OVERSAMPLING

               simcado.run(my\_src, ….. , SIM\_OVERSAMPLING=4)

If SimCADO is using the 0.004 arcsec mode, then it will calculate
everything based on a 0.001 arcsec grid. It resamples back up to 0.004
when passing the image to the detector module. Depending on the level of
accuracy you need (and the computer power you have available), you can
oversample as much as you’d like. Be careful with this - it uses lots of
RAM.

-  SimCADO also has a sub-pixel mode built into the method
   ``<Source>.apply_optical_train(... ,sub_pixel=False)``, - (very slow
   if there are >100 stars, it doesn’t use FFTs) - however I’ve only
   done a quick tests to make sure the code works. I can’t guarantee
   it’s accurate though. I’ll put that on my list of things to do this
   week (22 Nov 2016)

What SimCADO can do?
--------------------

What SimCADO can’t yet do?
--------------------------

What SimCADO will never do?
---------------------------

I have useful instrument data, who do I give it to?
---------------------------------------------------
