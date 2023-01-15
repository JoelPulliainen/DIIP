Test.m
Test estim_coins on all the measurement images.
Includes Imean which takes reads all the given calibration images and calculates means.
To use Imean all bias,dark and flatfield images need to be in the working directory and they all must have same names that were given out with the assignment.

estim_coins.m
Detect and Count coins in a image
Using Intensity calibration, geometrical calibration, removing checkerboard from the image, grayscaling, and filtering and filling holes of a image.
Input is Measurement image which contains coins to be counted
Mean bias image, mean dark image and flatfield image.
Output is s a six-element vector with the first element
corresponding to the number of 5 cent-coins in the image and the last element corresponding to
the number of 2 euro-coins in the image
Additional figures for visualising the process are included as comments.