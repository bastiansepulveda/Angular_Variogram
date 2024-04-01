functions.R

	File containing functions for working with images and estimators for angular variogram, cross angular variogram, and pseudo cross angular variogram.

	LIBRARIES
	The following libraries are called by the functions.R file

		fields, MASS, imager, writexl, progress, tidyverse

	FUNCTIONS
	Functions for processing images.
	- openImage(path):
		* parameters:
			path: the path of the image to be opened.
		* return:
			A matrix array that stores the grayscale pixel values.

	- rotateImage(img, ang):
		* parameters:
			img: square image to be rotated. Must be a square matrix array.
			ang: rotation angle.
		* return:
			A square matrix array with the rotated image, maintaining the original image's dimension. Corners of the image will be cropped if necessary.

	- scaleImage(img, l):
		* parameters:
			img: square image to be rescaled. Must be a square matrix array.
			l  : number of pixels on one side of the height (width) of the rescaled image.
		* return:
			A square matrix array with the rescaled image.

    Estimators.
    - VarAng(img, sep = 1):
        * parameters:
            img           : square image for which the angular variogram will be calculated. Must be a square matrix array.
            sep (optional): grid resolution. Must be an integer divisor of 180. If it is n, it calculates every n steps, from 1. Default is 1.
        * return:
            An array of size 180/sep with the angular variogram values at each considered angle.

    - VarAngCr(img1, img2, sep = 1):
        * parameters:
            img1, img2    : square images for which the cross angular variogram will be calculated. Must be square matrix arrays.
            sep (optional): grid resolution. Must be an integer divisor of 180. If it is n, it calculates every n steps, from 1. Default is 1.
        * return:
            An array of size 180/sep with the cross angular variogram values at each considered angle.

    - PseudoVarAngCr(img1, img2, sep = 1):
        * parameters:
            img1, img2    : square images for which the pseudo cross angular variogram will be calculated. Must be square matrix arrays.
            sep (optional): grid resolution. Must be an integer divisor of 180. If it is n, it calculates every n steps, from 1. Default is 1.
        * return:
            An array of size 180/sep with the pseudo cross angular variogram values at each considered angle.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

create_simulated_images.R

    This file is used to simulate images with a given covariance structure.

    LIBRARIES
    The following libraries are called by the create_simulated_images.R file

        fields, MASS, writexl, progress, rstudioapi

    FUNCTIONS
    Exponential, spherical, and Gaussian covariance functions. Each one has its parameters. It's possible to define new covariance functions.

    Function to simulate images.
    - simulateImages(nx, ny, N, covar, name_output_file):
        * parameters:
            nx              : number of pixels in the width of the simulated image.
            ny              : number of pixels in the height of the simulated image.
            N               : number of images to simulate.
            covar           : covariance function to use. Must be the name of a previously defined covariance function, in character format.
            name_output_file: name of the file that will contain the simulated images. Must be in character format.
        * return:
            A .csv file with N columns, each one corresponding to a simulated image.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

simulated_images.R

    This file calculates the cross angular variogram and pseudo cross angular variogram for the simulated images, delivering a summary of the results.

    LIBRARIES
    The following libraries are called by the simulated_images.R file

        rstudioapi

    FUNCTIONS
    - MeanStdResults(func, ang, images):
        * parameters:
            func   : functions to apply to the images. Must be a list with the names of the functions in character format. Can be "VarAngCr", "PseudoVarAngCr", or both.
            ang    : rotation angles of the image with respect to its center. Must be a list.
            images : images that will be used for the calculations. Must have the format of the output from the simulateImages function.
        * return:
            A matrix array with the function func used, the angle, the mean of the results, and the standard deviations.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

real_images.R

    This file calculates the cross angular variogram and pseudo cross angular variogram of real images, delivering a summary of the results.

    LIBRARIES
    The following libraries are called by the real_images.R file

        rstudioapi
