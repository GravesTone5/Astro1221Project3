# Astro1221Project3

Instructions:
Fit star shapes in images to measure brightnesses, then create a color-magnitude diagram to study stellar evolution.

Astronomy Context: Open clusters are groups of stars born together from the same gas cloud, so they're all the same age and distance. A color-magnitude diagram (CMD) plots each star's color (related to temperature) vs. brightness. The pattern reveals stellar evolution: hot blue stars on upper-left, Sun-like stars in the middle "main sequence," red giants above. By comparing to theoretical models, we determine cluster age and distance. PSF (Point Spread Function) fitting means fitting each star's image to extract its brightness precisely.

What You'll Do: Download HST or ground-based images of an open cluster (M67, NGC 2264) from MAST. Use Astropy to read the FITS file and display the image. Identify star positions (you can do this manually by clicking on stars, or use simple source detection). For each star, extract a small cutout (like 15×15 pixels). Fit a 2D Gaussian: I(x,y) = A exp[-((x-x₀)²+(y-y₀)²)/(2σ²)] + background. The amplitude A tells you the star's brightness.

Technical Note on PSF Fitting: If you fix σ (the PSF width) from bright, well-measured stars, then flux ∝ A directly—this is the simple approach. If you let σ vary per star, you need to calculate flux = 2πAσ². For this project, fixing σ is recommended! After fitting, convert flux to magnitude using m = -2.5 log₁₀(flux) + constant.

Color-Magnitude Diagram: If you have images in two filters (e.g., B and V, or g and r), fit stars in both images. The color is the magnitude difference (B-V or g-r). Plot color vs. magnitude to make a CMD. If you only have one filter, you can still make a magnitude-vs-position plot or single-filter CMD.

About Magnitudes: You'll get "instrumental magnitudes" (relative brightnesses) unless you calibrate with standard stars. That's fine! The CMD shape is what matters. For distance estimates, you'd need to either calibrate your magnitudes or compare to theoretical models (isochrones) that you shift vertically to match your CMD.

Success Looks Like: Fitting 10-15 stars in an image, implementing 2D Gaussian PSF model with fixed σ, using scipy.optimize.curve_fit or least_squares for each star, converting to magnitudes, creating a CMD (or magnitude distribution), and discussing what you see (main sequence? red giants? cluster age estimate?).

Advanced Options: Fit multiple filters for true color-magnitude diagram. Implement background subtraction. Handle crowded fields with overlapping PSFs. Measure σ from bright stars. Apply aperture corrections. Compare to theoretical isochrones to determine age and distance.