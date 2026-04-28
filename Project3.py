import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from astropy.io import fits
import scipy 
import urllib.request
import os 
from scipy import optimize
from matplotlib.colors import LogNorm
import array
import os 
from scipy import optimize
import photutils.detection as detect

NGC2660_555_fits = "hst_10634_03_acs_wfc_f555w_j9dm03_drc.fits"
NGC2660_814_fits = "hst_10634_03_acs_wfc_f814w_j9dm03_drc.fits"
NGC2660_fits = "hst_10634_03_acs_wfc_total_j9dm03.fits"
NGC2660_fits_url = "https://hla.stsci.edu/cgi-bin/getdata.cgi?config=ops&dataset=hst_10634_03_acs_wfc_total_j9dm03"
NGC2660_555_fits_url = "https://hla.stsci.edu/cgi-bin/getdata.cgi?config=ops&download=1&dataset=hst_10634_03_acs_wfc_total_j9dm03&filename=hst_10634_03_acs_wfc_f555w_j9dm03_drc.fits"
NGC2660_814_fits_url = "https://hla.stsci.edu/cgi-bin/getdata.cgi?config=ops&download=1&dataset=hst_10634_03_acs_wfc_total_j9dm03&filename=hst_10634_03_acs_wfc_f814w_j9dm03_drc.fits"

if not os.path.exists(NGC2660_814_fits) & os.path.exists(NGC2660_555_fits):
        try:
            urllib.request.urlretrieve(NGC2660_814_fits_url, NGC2660_814_fits)
            print(f"Downloaded {NGC2660_814_fits}")
            urllib.request.urlretrieve(NGC2660_555_fits_url, NGC2660_555_fits)
            print(f"Downloaded {NGC2660_555_fits}")
        except Exception as e:
            print("Error Downloading File")
            NGC2660_555_fits = None
            NGC2660_814_fits = None
bands = [NGC2660_555_fits, NGC2660_814_fits]
image_data = []
for b in bands:
    with fits.open(b) as HDUList:
        image_data.append(HDUList[1].data) #reads the data from the SCI section of the fits file
        HDUList.close()
image_555, image_814 = image_data

sigma = 3
fwhm = 2.355 * sigma
clipped_data_555 = image_555.copy()
for i in range(5):
    mean = np.nanmean(clipped_data_555)
    std = np.nanstd(clipped_data_555)
    mask = np.abs(clipped_data_555 - mean) <= 3 * std
    clipped_data = clipped_data_555[mask]
    background_mean_555 = np.mean(clipped_data)

clipped_data_814 = image_814.copy()
for i in range(5):
    mean = np.nanmean(clipped_data_814)
    std = np.nanstd(clipped_data_814)
    mask = np.abs(clipped_data_814 - mean) <= 3 * std
    clipped_data = clipped_data_814[mask]
    background_mean_814 = np.mean(clipped_data)

subtracted_image_555 = image_555 - background_mean_555

plt.figure()
plt.imshow(subtracted_image_555, cmap='grey', origin = 'lower', norm=LogNorm())
plt.colorbar(norm=LogNorm())
plt.show()

fig, axes = plt.subplots()
axes.hist(image_555.flatten(), bins=50, color='blue', range=(-.1, .75), alpha=.6)
axes.hist(subtracted_image_555.flatten(), bins=50, color='red', range=(-.1, .75), alpha=.6)

subtracted_image_814 = image_814 - background_mean_814

plt.figure()
plt.imshow(subtracted_image_814, cmap='grey', origin = 'lower', norm=LogNorm())
plt.colorbar(norm=LogNorm())
plt.show()

fig, axes = plt.subplots()
axes.hist(image_814.flatten(), bins=50, color='blue', range=(-.1, .75), alpha=.6)
axes.hist(subtracted_image_814.flatten(), bins=50, color='red', range=(-.1, .75), alpha=.6)

subtracted_image_555 = image_555 - background_mean_555
subtracted_image_814 = image_814 - background_mean_814
detection_threshold = 5 * background_mean_555

filter = int(np.ceil(fwhm))

local_max = scipy.ndimage.maximum_filter(subtracted_image_555, size=filter)

peaks = (subtracted_image_555 == local_max) & (subtracted_image_555 > detection_threshold)

detections_y, detections_x = np.where(peaks)

fig, axes = plt.subplots()

plt.figure()
plt.imshow(subtracted_image_555, cmap='grey', origin = 'lower', norm=LogNorm(vmin = .0001, vmax=300))
axes.scatter(detections_x, detections_y, color = 'purple', alpha =.6, marker='x')
plt.colorbar(norm=LogNorm(vmin = .0001, vmax=300))
plt.show()

class Star:
    """One star is the coordinates of the center of a star on the image"""
    
    def __init__(self, row, col, name=""):
        self.row = row
        self.col = col
        self.name = name

    def get_row_col(self):
        return (self.row, self.col)

    def value_at(self, image):
        return image[self.row, self.col]

    def cutout(self, image, half_size=8):
        nrows, ncols = image.shape
        r_start = self.row - half_size
        r_end = self.row + half_size + 1
        c_start = self.col - half_size + 1
        c_end = self.col + half_size + 1
        if r_start < 0:
            r_start = 0
        if c_start < 0: 
            c_start = 0
        if r_end > nrows:
            r_end = nrows
        if c_end > ncols:
            c_end = ncols
        return image[r_start:r_end, c_start:c_end]


stars = [
    Star(1311, 1231, "1"),
    Star(544, 3297, "2"),
    Star(2910, 4128, "3"),
    Star(4462, 1557, "4"),
    Star(1657, 2605, "5"),
    Star(4064, 1437, "6"),
    Star(2160, 690, "7"),
    Star(4531, 1714, "8"),
    Star(2010, 3225, "9"),
    Star(2710, 1710, "10"),
    Star(4531, 1713, "11"),
    Star(4168, 2647, "12"),
    Star(3565, 3130, "13"),
    Star(3203, 1623, "14")
]


flux = []

def psf_fit(coords, row, col, sigma=3, amplitude_func=100, offset=0):
    x, y = coords
    

    psf = amplitude_func * np.exp(-((y-row)**2 + (x-col)**2)/(2*sigma**2))
    fwhm = 2.355*sigma
    return psf.ravel()

def image_2660(selected_image):
    for s in stars:
        cutout = s.cutout(image=selected_image)
        y_size, x_size = cutout.shape
        yy, xx = np.meshgrid(np.arange(y_size), np.arange(x_size), indexing='ij') 
        xdata = (yy, xx)
        ydata = cutout.ravel()
        print(ydata.shape)
        ydata = np.nan_to_num(ydata, nan=np.nanmedian(ydata))
        peak = np.max(ydata)
        bg = np.median(ydata)
        flux.append(2 * np.pi * peak * 2**2)

        lower_bounds = [0, 0, 0.5, 0, -np.inf]
        upper_bounds = [y_size, x_size, 10.0, np.inf, np.inf]
        
        p0 = [
        y_size / 2,
        x_size / 2,
        1.0,
        peak - bg,
        bg
        ]
        print(p0)
        popt, pcov = optimize.curve_fit(
            psf_fit,
            xdata,
            ydata,
            p0=p0,
            bounds=(lower_bounds, upper_bounds),
            maxfev=2000
        )
        fit_model = psf_fit(xdata, *popt).reshape(y_size, x_size)
        actual_data = ydata.reshape(y_size, x_size)
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))

        im0 = ax[0].imshow(actual_data, origin='lower')
        ax[0].set_title('Actual Star (Data)')
        fig.colorbar(im0, ax=ax[0])
        
        im1 = ax[1].imshow(fit_model, origin='lower')
        ax[1].set_title('Best Fit')
        fig.colorbar(im1, ax=ax[1])
        
        im2 = ax[2].imshow(actual_data - fit_model, origin='lower')
        ax[2].set_title('Residuals (Data - Model)')
        fig.colorbar(im2, ax=ax[2]) 
        print(f"popt {popt}")
image_2660(subtracted_image_555)
image_2660(subtracted_image_814)
print(f"flux {flux}")

flux_split = len(flux)//2 
intensity_555 = np.array(flux[:flux_split])
intensity_814 = np.array(flux[flux_split:])
print(f"flux_555 {intensity_555}")
print(f"flux_814 {intensity_814}")
flux_subtraction = intensity_555 - intensity_814
print(f"flux_subtraction {flux_subtraction}")

plt.scatter(flux_subtraction, intensity_814)
plt.xlabel("Color/Temperature")
plt.ylabel("Magnitude/Luminosity")
plt.title("NGC 2660 Color-Magnitude Diagram")
plt.show()