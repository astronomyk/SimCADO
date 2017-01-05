# Updates to SimCADO

## Data updates
The data files used by SimCADO are continually being updated to reflect the newest information coming from the MICADO work packages. 

To update your version of SimCADO, use the `simcado.get_extras()` command.

### 2017-01-05
* VLT mirror coating
        TC_aluminium.dat        
* HAWK-I filter curves: 
        TC_filter_J.dat         
        TC_filter_H.dat         
        TC_filter_Ks.dat        
        TC_filter_Y.dat         
* Updated detector layout
        FPA_chip_layout.dat     


### 2016-11-12
* `TC_surface.dat` 
  Updated with more optimistic Strehl values for JHK central wavelengths
* `TC_mirror_gold.dat`
  Added a reflectivity curve for an unprotected gold coating
* `default.config`
  Changed `INST_MIRROR_TC` to reference `TC_mirror_gold.dat`


## Source Code updates
The current stable version is SimCADO v0.3. The current development version is SimCADO v0.4dev

### Development version
2017-01-05
* SimCADO now preloads the transmission curves. Generating an OpticalTrain now only takes ~2 seconds (compared to 20 previously)

2016-12-06
* Added functionality to generate "ideal" AO PSFs using POPPY for the diffraction limited core and added a Seeing halo

2016-11-12
* Bug fix in simcado.psf - psf.lam_bin_centers was taking opt_train.lam_bin_centers instead of using the "WAVE0" keywords in the FITS header
* Bug fix in source.psf - apply_optical_train() was asking for PSFs outisde of the opt.train.psf range. related to the point above