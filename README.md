# Astro1221 Project 3

> ## Siena Baez, Xander Graves, Alexis Walker

# Topic 8: Open Cluster Photometry and Color-Magnitude Diagram
> ## Goal:
>> The goal is to take an open cluster (NGC 2660 in our case) from MAST and collect fits files in order to gather information about the stars in the cluster. After that, we needed to find the pixel intensity (of both bands) as well as the centers of the stars in order to use PSF to get the colors and the magnitudes. The data was taken from HST (Hubble Space Telescope) and then put into a URL to store it, and then we used Astropy to open the file and get the information of the cluster.

> ## Requirements and Use of Packages:
>> Looking into Requirements.txt, you will see that you need to install 2 or 3 different packages (Astropy, Photutils and Scipy). All of these are for extracting or gathering information from the HDU of the fits file. 

>> Astropy: For this one, you will see that you need it for openining the fits file in the first place. You will need files of both the bands, in order to get the PSF fittings and the color-magntiude diagram. This is optional (though I find it cool) but you can also open the total (the original image given) fits file to see the original image with both the bands on it.

>> Photutils: This is optional, as we were never able to actually implement it our code. We wanted to use it for source detections (or finding out where stars are), but it was not working the way we wanted it to. We are sure there are ways that work better than our method, so feel free to try!

>> Scipy: If not already included in your python distribution, !pip install it in order to use it for optimization for the curve-fitting of the stars.

> ## Key Takeaways:
>> We learned that with our image that HST used different "bands" to act as filters to get the blue and orange colors of the stars. The catch to this was that they did not use blue or orange wavelengths for the image, they used something called "false imaging." THis meant that the blue was actually more of a yellow-green (555 nm) and the orange was more of a infrared wavelength (814 nm). These were the two bands that we used to check for the intensity and the colors of the stars (which we will get into more in the docstrings).

>> On top of this, we also learned more about how color-magnitude diagrams work, and how to better extract information from the bands. With what we were given, one of the wavelengths (555 nm) was used to be the intensity of the stars (or the y-axis on the CMD) and for the x-axis, we subtracted the wavelength we did not use (814 nm) from the one we did use.

> # AI Usage on This Project
>> We used AI to help us with a lot of the issues we had with the PSF fitting (especially if we could not meet with any of the instructors). This was mainly to correct errors as well as point us into the right direction to finish up our code. We also used it to help with the Star class (as none of us had any prior knowledge on how to make classes). It also helped near the beginning when we were starting to open the fits files using Astropy. We tried to keep AI usage as minimal as possible; however, with all the errors we had no idea how to fix, we used it as a guide. When we did not use AI, we researched with either the package's website or watching videos of people who know a lot more about code than we do.