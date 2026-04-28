import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astropy.units as u
import astropy.constants as const
from astropy.io import fits
import scipy 
import urllib.request
import os 
from scipy import optimize
import photutils.detection as detect
import streamlit as st

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

HDUList = fits.open(NGC2660_fits)

image_total = HDUList[1].data

plt.figure()
plt.imshow(image_total, cmap='gray', origin = 'lower', norm=LogNorm())
plt.colorbar(norm=LogNorm())
plt.show()

plt.savefig('NGC2660_total.png')

HDUList.close()

bands = [NGC2660_555_fits, NGC2660_814_fits]
image_data = []
for b in bands:
    with fits.open(b) as HDUList:
        image_data.append(HDUList[1].data) 
        HDUList.close()
image_555, image_814 = image_data

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

subtracted_image_555 = image_555 - background_mean_555

plt.figure()
plt.imshow(subtracted_image_555, cmap='grey', origin = 'lower', norm=LogNorm())
plt.colorbar(norm=LogNorm())
plt.show()

fig, axes = plt.subplots()
axes.hist(image_555.flatten(), bins=50, color='blue', range=(-.1, .75), alpha=.6)
axes.hist(subtracted_image_555.flatten(), bins=50, color='red', range=(-.1, .75), alpha=.6)

sigma = 3

fwhm = 2.355 * sigma

subtracted_image_555 = image_555 - background_mean_555
subtracted_image_814 = image_814 - background_mean_814
print(background_mean_555)
detection_threshold =  5 * background_mean_555
print(detection_threshold)
filter = int(np.ceil(fwhm))


detections =  detect.find_peaks(data=subtracted_image_555, threshold=detection_threshold, mask=mask*1000, min_separation=50 , n_peaks=200)
print(detections)
detections_x = detections['x_peak']
detections_y = detections['y_peak']

plt.figure()
plt.imshow(subtracted_image_555, cmap='grey', origin = 'lower', norm=LogNorm(vmin = .0001, vmax=300))
plt.scatter(detections_x, detections_y, color = 'purple', alpha =.6, marker='o', )
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
    Star(2010, 3230, "9"),
    Star(2710, 1710, "10"),
    Star(4531, 1713, "11"),
    Star(4168, 2647, "12"),
    Star(3565, 3130, "13"),
    Star(3203, 1623, "14"),
    Star(2228, 2120, "15"),
    Star(2888, 2795, "16"),
    Star(2015, 2177, "18"),
    Star(3056, 1063, "19"),
    Star(4116, 3859, "20"),
    Star(3333, 866, "21"),
    Star(2340, 1630, "22"),
    Star(1880, 2087, "23"),
    Star(1900, 1738, "24"),
    Star(2003, 1967, "25"),
    Star(1953, 1999, "26"),
    Star(1988, 2106, "27"),
    Star(2016, 2176, "28"),
    Star(2193, 2216, "29"),
    Star(1433, 2298, "30")
    ]

flux = []

def psf_fit(coords, row, col, sigma=3, amplitude_func=10, offset=0):
    x, y = coords

    psf = amplitude_func * np.exp(-((y-row)**2 + (x-col)**2)/(2*sigma**2))
    fwhm = 2.355*sigma
    return psf.ravel()

def image_2660(selected_image):
    for s in stars:
        cutout = s.cutout(image=selected_image)
        nonan_cutout = np.nan_to_num(cutout)
        y_size, x_size = nonan_cutout.shape
        yy, xx = np.meshgrid(np.arange(y_size), np.arange(x_size), indexing='ij') 
        coordinates = (yy, xx)
        amplitude = nonan_cutout.ravel()
        amplitude = np.nan_to_num(amplitude, nan=np.nanmedian(amplitude))
        peak = np.max(amplitude)
        bg = np.median(amplitude)
    

        lower_bounds = [0, 0, 0.5, 0, -np.inf]
        upper_bounds = [y_size, x_size, 10.0, np.inf, np.inf]
        
        p0 = [
        y_size / 2,
        x_size / 2,
        1.0,
        peak - bg,
        bg

        ]
        popt, pcov = optimize.curve_fit(
            psf_fit,
            coordinates,
            amplitude,
            p0=p0,
            bounds=(lower_bounds, upper_bounds),
            maxfev=20000
            )         
        actual_data = amplitude.reshape(y_size, x_size)
        first_fit_model = psf_fit(coordinates, *popt).reshape(y_size, x_size)
        mask_data = (actual_data - first_fit_model) > 60.
        actual_data[mask_data] = np.nan
        nonan_masked = np.nan_to_num(actual_data, nan=background_mean_555)
        masked_popt, mask_pcov = optimize.curve_fit(
            psf_fit,
            coordinates,
            nonan_masked.ravel(),
            p0=p0,
            bounds = (lower_bounds, upper_bounds),
            maxfev=500000,
            )

        actual_data = amplitude.reshape(y_size, x_size)
        peak_final = np.max(nonan_masked) 
        flux.append(2 * np.pi * peak_final * 2**2)
        
        
image_2660(subtracted_image_555)
image_2660(subtracted_image_814)
print(f"flux {flux}")

flux_split = len(flux)//2 
intensity_555 = np.asarray(flux[:flux_split])
intensity_814 = np.asarray(flux[flux_split:])
print(f"flux_555 {intensity_555}")
print(f"flux_814 {intensity_814}")
flux_subtraction = intensity_555 - intensity_814
print(f"flux_subtraction {flux_subtraction}")

plt.figure()
flux_plot = plt.scatter(flux_subtraction, intensity_814)
plt.xlim(-400, 200)
plt.ylim(5000, 6150)
plt.xlabel("Color")
plt.ylabel("Magnitude")
plt.title("NGC 2660 Color-Magnitude Diagram")
plt.show()

plt.savefig('NGC2660_CMD.png')

st.title("NGC 2660 Color-Magnitude Diagram")
st.write("The image of the total cluster is shown below:")
st.image('NGC2660_total.png')
st.write("This is the Color-Magnitude Diagram of NGC 2660, which is an open cluster of stars.")
st.write("Press the button below to see the color-magnitude diagram of manually selected stars:")
if st.button("See Color-Magnitude Diagram"):
    with st.spinner("Loading color-magnitude diagram..."):
        st.subheader("Color-Magnitude Diagram")
        st.image('NGC2660_CMD.png')