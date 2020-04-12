Sample images have been included for reference. 
Example Command
A2b_Paras_2018CSB1111_2020_CS517('13954.jpg','donald.tif',tpts,'test',1);

This would show the intermediate images and also output the result gif file as test.gif.

Approach:
I have basically used delaunay triangulation on tie points of the intermediate images which are generated using ratio*im1_pts+(1-ratio)*im2_pts
The using this delaunay triangulation I find out the corresponding triangle in im1 and im2. For all points within an triangle in the intermediate image 
I compute the barycentric coordinates with respect to the triangle containing it. Using the barcentric coordinates I try to find out the correspoinding points by coverting barycentric back to cartesian w.r.t the triangle coordinates in the im1 and im2.
When I find the coordinates I used the bilinear Interpolation to find out the pixel value in intermediate image from im1 and im2 and use the weighted value as im1*ratio +(1-ratio)*im2;
Thus by varying the ratio I computed all the intermediate images and put htem in a GIF file.


To generate the tie points I have used the dlib library in python which had a model for '68 facial landmarks' and after obatining the 68 facial tie points I dumped them to mat file which contains 68*4 tpts array to be used in the input to the matlab function. This ipyNoteook could be used to used to geenrate tie points for any two images and is best is ran on Google Colab as it conatins all the necessary libraries 
to generate the tie points.
The mat file could be loaded easily on matlab.

References:

http://www.seas.upenn.edu/~cse399b/Lectures/CSE399b-07-triangle.pdf
https://github.com/owenqyzhang/Face-Morphing
