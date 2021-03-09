Nanoparticle Sizing and Morphology Analysis Program
***********************************************
by Alex Powers  alex.powers@berkeley.edu 

v.2.0
1/26/15
2/10/15 changed window appearance location 


Introduction 
***********
This program works best with images containing spherical particles with large inter-particle spacing and high contrast. 
However, the algorithm is very capable of dealing with high signal to noise ratio, overlapping and touching particles, and low or uneven contrast. 
The algorithms are not optimized for rods and other elongated particles due to the implementation of the touching particle separation methods. 
For a set of sample .dm3 images for which this program performs within 3% deviation ( of average area) of outlining by hand, see the 'Test Images' folder. 
It is encouraged to try the following directions using those images to see how the program works on a reliable data set. 

Image Requirements 
******************
The program is currently compatible with .png .tif and .dm3 file types. 
It is important that the algorithm knows the pixel to nanometer ratio for each image. 
For .dm3 files, this information is extracted automatically. 
For .tif or .png files the scale needs to be appended to the image file name after a '$' 
Example: if for 'sample1.tif' there are 10.0 pixels / nanometer, rename the file 'sample1$10.0.tif'
The pixel / nm scale can be calculated by determining the pixel length of a scale bar in a program like imagej or photoshop.  

A set of images corresponding to a particular sample should be grouped together in a folder for analysis. 
Different samples should be in separate folders. One sample or folder can be analyzed at a time. 
All images in the folder you wish to analyze should be of the same filetype. 
To analyze another sample; close all the figures and re run the program. 

Instructions
**********
To start the program, double click  the gui_analysis.m file to open in MATLAB. 
Ensure that the folder containing this file is now the current folder or in the MATLAB path. 
Click Run for the gui_analysis.m file to start the interactive user interface. 

3 new figure windows should open. The left most figure contains the 
control panel and the other two display results. 

If you wish to use images with a black background, check the appropriate box in the control panel. 

Now, click 'select folder' button and in the pop up window choose the folder corresponding to your images and select it. 
The analysis should proceed automatically. 
The progress is displayed on the MATLAB command line - the algorithm takes ~8 seconds per 2000 x 2000 pixel image. 
If any errors result, they may be displayed in the command line. Please contact me. 

After the initial analysis is complete, the analyzed images should appear in figure 3. The particles found are outlined. 
Because large clumps of particles are also selected, further filtering is necessary to select only separate individual particles. 
These filters are provided in the control panel. The automatic filtering is to cutoff particles with circularity (also known as sphericity) below 0.7. 
For a image of sphericity vs roundness see http://www.globalenergylaboratories.com/wp-content/uploads/2013/08/roundness-sphericity-matrix.jpg. 

'true' particles that meet the filter criteria are outlined in red. 'false' particles are outlined in blue. 
Below each image in figure 3, the threshold used to the select the particles in a particular image is displayed; this is a number x from 0 - 1 of which all 
pixels with a brightness below x are selected. A threshold algorithm is used that was optimal for a set of test images - this may not be the best 
in each individual case due to differences in image contrast and brightness. 

A new threshold can be applied (just for that image); enter the new threshold and click enter - the image and results will update. 
In figure 3, there is a drop down menu from which different  images can be selected for display and modification. 

For an example of interactively changing the threshold, try analyzing the images in the test images folder. 
For Image 1, the predicted threshold is 0.3 (due to low contrast for most particles) whereas the threshold for greatest accuracy is ~ 0.6.
 Note that the large clumps are filtered out by circularity. 
In contrast (pun intended) image 2 in the test images folder overestimates the threshold at 0.6 where actually 0.45 is optimal. 
The low magnification in image 2 makes it non optimal for high accuracy sizing. 
The threshold determined for Image 3 and 4 delivers results within 5% error of hand selecting and only slight adjustment is necessary if high accuracy is required. 

Changing the filtering immediately updates the 2 histograms, averages, and SDs of size in figure 1. The third plot in figure 1 contains 
a spread of the ALL particles (true and false) by roundness and circularity with coloring based on diameter (from area). 

Provide the sample name in the box of the control panel to add this to the histograms. 
Clicking save fig. or save .tiff will save figure 1 as a file in the directory containing the images. 



