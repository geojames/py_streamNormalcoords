# py_streamNormalcoords

Stream Normal Coordinate Transformation
Compatibility:  Python 3.7

This program is based on an academic research article:
> Legleiter CJ, Kyriakidis PC. 2006. Forward and Inverse Transformation  between Cartesian and Channel-fitted Coordinate Systems for Meandering Rivers. Mathematical Geology 38 : 927–958. DOI: 10.1007/s11004-006-9056-6
>- The original source code from this article was in Matlab                 (acquired by Dietrich from Carl Legleiter, cjl@usgs.gov). This code is a translation of that original Matlab code with some added features/tweaks.
                
Given a digitized stream centerline and a set of input x,y         coordinates. This software will transform the input x,y coordinates into a curvalinear coordinate system based on the stream centerline. This stream normal coordinate system is:
- **ds** - a downstream distance (from the most upstream point)
- **xs** - a cross-stream distance (an orthogonal distance to the centerline along the normal vector)
