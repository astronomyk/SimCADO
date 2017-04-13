import numpy as np
import matplotlib.pyplot as plt

from photutils import Background2D, SigmaClip, MedianBackground
from photutils import DAOStarFinder
#from photutils import aperture_photometry, CircularAperture, CircularAnnulus
from astropy.stats import sigma_clipped_stats

class Stamp(object):

    def __init__(self, image, image_id=None, **kwargs):

        self.params = {"id"     : image_id,
                       "mag_zero_point"    : 0,
                       "aperture_radius"   : 3.,
                       "sky_radius_inner"  : None,
                       "sky_radius_outer"  : None,
                       "hot_pixel_threshold" : 3000
                      }
        self.params.update(kwargs)

        self.image = remove_hot_pixels(image, threshold=self.params["hot_pixel_threshold"])
        
        #try:
        #   self.bg_mean, self.bg_median, self.bg_std = \
        #                    sigma_clipped_stats(self.image, sigma=3., iters=3)
        #except:
        #    self.bg_mean, self.bg_median, self.bg_std = 0, 0, 0

        self.mag, self.snr, self.flux, self.noise, self.bg_std, self.bg_median = \
                            self.get_photometry(self.params["aperture_radius"],
                                                self.params["sky_radius_inner"],
                                                self.params["sky_radius_outer"],
                                                self.params["mag_zero_point"])

        try:
            height, x, y, width_x, width_y = self.get_moments()
        except:
            height, x, y, width_x, width_y = [0]*5

        self.peak = height
        self.x,  self.y  = x, y
        self.dx, self.dy = width_x, width_y
        self.fwhm = self.get_fwhm()

    def get_photometry(self, aperture_radius=3,
                       sky_radius_inner=None, sky_radius_outer=None,
                       mag_zero_point=0, mode="basic"):


        if sky_radius_outer is None:
            sky_radius_outer = np.min(self.image.shape) // 2
        if sky_radius_inner is None:
            sky_radius_inner = sky_radius_outer - 3

        x, y = np.array(self.image.shape) // 2

        # use photutils?
        # aperture = CircularAperture((x, y), r=aperture_radius)
        # bg_aperture = CircularAnnulus((x, y), r_in=sky_radius_inner,
        #                                      r_out=sky_radius_outer)
        # phot_table = aperture_photometry(image, (aperture, bg_aperture))

        # bg = phot_table["aperture_sum_1"] / bg_aperture.area()
        # flux = phot_table["aperture_sum_0"] - bg * aperture.area()

        # noise = self.std

        r  = aperture_radius
        ro = sky_radius_outer
        dw = sky_radius_outer - sky_radius_inner

        if "circ" in mode.lower():
            pass
            
            
        elif "basic" in mode.lower():

            bg = np.copy(self.image[y-ro:y+ro, x-ro:x+ro])
            bg[dw:-dw, dw:-dw] = 0
            bg_median = np.median(bg[bg!=0])
            bg_std = np.std(bg[bg!=0])

            im = np.copy(self.image[y-r:y+r, x-r:x+r])
            flux = np.sum(im - bg_median)

            noise = bg_std * np.sqrt(np.sum(bg != 0))
            snr = flux / noise

        mag = -2.5 * np.log10(flux) + mag_zero_point

        return mag, snr, flux, noise, bg_std, bg_median


    def get_moments(self, radius=None):
        
        #x_cen, y_cen = self.image.shape // 2
        
        #if np.isscalar(radius):
    
        return gaussian_moments(self.image)

    def get_fwhm(self):
        
        return 2.35 * np.sqrt(self.dx**2 + self.dy**2)


class PostageStamps(object):

    def __init__(self, image, x=None, y=None, name=None, **kwargs):


        self.params = {"bg_tile_size"      : 32,
                       "fwhm"              : 5.,
                       "stamp_width"       : 32,
                       "stamp_height"      : None,
                       "mag_zero_point"    : 0,
                       "aperture_radius"   : 3.,
                       "sky_radius_inner"  : None,
                       "sky_radius_outer"  : None,
                       "threshold"         : None,
                      }

        self.params.update(kwargs)


        self.image_bg = self._get_background(image,
                                             tile_size=self.params["bg_tile_size"])
        self.image = image - self.image_bg

        if x is None and y is None:
            self.x, self.y = self.find_sources(fwhm=self.params["fwhm"],
                                               threshold=self.params["threshold"])
        elif x is not None and y is not None:
            self.x, self.y = x, y
        else:
            raise ValueError("x and y need to be both None or equal length arrays")

        self.stamps = self.get_stamps(self.x, self.y,
                                      w=self.params["stamp_width"],
                                      h=self.params["stamp_height"])

        self._get_photometry()


    def find_sources(self, fwhm=5., threshold=None):

        mean, median, std = sigma_clipped_stats(self.image, sigma=3., iters=3)

        if threshold is None:
            threshold=5.*std

        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
        sources = daofind(self.image - median)
        self.DAOsources = sources

        return sources["xcentroid"], sources["ycentroid"]


    def get_stamps(self, x, y, w=16, h=None):

        if h is None:
            h = w

        x0, x1 = x - w//2, x - w//2 + w
        y0, y1 = y - h//2, y - h//2 + h
        ims = [self.image[yy0:yy1, xx0:xx1] for yy0, yy1, xx0, xx1 in zip(y0, y1, x0, x1)]

        stamps = [Stamp(im, **self.params) for im in ims]

        return stamps


    def _get_background(self, image, tile_size=32):

        sigma_clip = SigmaClip(sigma=3., iters=3)
        bkg_estimator = MedianBackground()
        bkg = Background2D(image, tile_size, filter_size=3,
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        return bkg.background


    def _get_photometry(self):

        self.mags      = np.array([stamp.mag   for stamp in self.stamps])
        self.snrs      = np.array([stamp.snr   for stamp in self.stamps])
        self.fluxes    = np.array([stamp.flux  for stamp in self.stamps])
        self.noises    = np.array([stamp.noise for stamp in self.stamps])
        self.fwhms     = np.array([stamp.fwhm  for stamp in self.stamps])
        self.peaks     = np.array([stamp.peak  for stamp in self.stamps])
        self.bg_stds   = np.array([stamp.bg_std     for stamp in self.stamps])
        self.bg_medians = np.array([stamp.bg_median for stamp in self.stamps])
        

    def plot_stamps(self, n, n_wide=5, colorbar=False, vmin=None, vmax=None, norm=None):

        if isinstance(n, str) and n == "all":
            n = range(len(self.stamps))

        if np.isscalar(n):
            n = range(n)

        w = n_wide
        l = len(n)
        h = l // w + 1

        for i in n:
            plt.subplot(h, w, i+1)
            plt.imshow(self.stamps[i].image, norm=norm, vmin=vmin, vmax=vmax)

            if colorbar:
                plt.colorbar()


def gaussian_moments(data):
    """
    Returns (height, x, y, width_x, width_y) for a 2D Gaussian

    Taken from:
        http://scipy-cookbook.readthedocs.io/items/FittingData.html

    Parameters
    ----------
    data : np.array
        2D Gaussian

    Returns
    -------
    height, x, y, width_x, width_y

    """
    height = data.max()

    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total

    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())

    return height, x, y, width_x, width_y
    
    
def remove_hot_pixels(image, threshold=5000):    
    """
    
    Parameters
    ----------
    image : np.ndarray
        The image, in 2D format
    threshold : float, int
        The difference between any two neighbbouring pixels. 

    
    Returns
    -------
    a : np.array
        The input image with the hot pixels set to the image median
    
    """
    a = np.copy(image)
    b = a[:-1,:] - a[1:,:]
    c = b[1:,:] - b[:-1,:]

    x, y = np.where(c > threshold)
    a[x+1, y] = np.median(a)
    
    return a