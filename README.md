
# py_streamNormalcoords

Stream Normal Coordinate Transformation
Compatibility:  Python 3.7

This program is based on an academic research article:
> Legleiter CJ, Kyriakidis PC. 2006. Forward and Inverse Transformation  between Cartesian and Channel-fitted Coordinate Systems for Meandering Rivers. Mathematical Geology 38 : 927–958. DOI: 10.1007/s11004-006-9056-6
>- The original source code from this article was in Matlab                 (acquired by Dietrich from Carl Legleiter, cjl@usgs.gov). This code is a translation of that original Matlab code with some added features/tweaks.
                
Given a digitized stream centerline and a set of input x,y         coordinates. This software will transform the input x,y coordinates into a curvalinear coordinate system based on the stream centerline. This stream normal coordinate system is:
- **ds** - a downstream distance (from the most upstream point)
- **xs** - a cross-stream distance (an orthogonal distance to the centerline along the normal vector)

![enter image description here](https://imgur.com/2MyvFN4.png)
        
### Data Prep
All input files are CSV format
##### Centerline Point (CL_pts)
X,Y points describing the centerline 
Minimum columns names = X, Y (caps)
> a spatially variable search radii can be added by adding extra columns to the input centerline data
>>- downstream variable radii, but symmetric - add 'R' column with distances
>>- downstream variable, but non-symmetric - add 'LR' and 'RR' for different left and right distances
    
##### Data Point (data_pts)
X,Y data points to be transformed
Minimum columns = X, Y (caps)
- additional columns will be transferred to the output

### To Run:
Ensure you have the following packages installed:
- scipy, numpy, pandas, matplotlib

Open the xy2sn.py file in a Python editor
Fill in the input and output file paths/names in the Input section
#### Choose your Transformation Parameters
The transformation parameters are important, but require some trial and error. A lot will depend on what your initial point spacing is on your centerline and how tight your meander bends are.
- nFilt  = number of filtering iterations
- order  = polynomial order of the filter
- window = number of points to include in the filter window
- nDiscr = number of segments to split the centerline into
	- a good place to start is your total length / desired segment length
- rMax   = maximum distance away from the centerline the code will search for points (see above how to make this spatially variable)

Run the code (use the run button in your editor)

### More Theory (coming soon)
- More in-depth description of the process
- Overlapping search radii in sharp corners
	- Spatially variable rMax values helps
	- With Uniform RMax there is the possibility of overlapping cells
	- ![enter image description here](https://imgur.com/DXKeiwa)
	- With Variable RMax there overlapping cells can be avoided
	- ![enter image description here](https://imgur.com/SNqJkAY)
- Downstream distance variability as a function of transform parameters
