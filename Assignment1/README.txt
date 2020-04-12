A2a_Paras_2018CSB1111_2020_CS517(image_name,qID,angle,dummy2,afnTM,output,toshow)

image_name: Name of the input Image along with an extension 

qID:	1 if angle needs to be specified
	2 if afnTM(affine matrix) needs to be specified

dummy2: dummy Variable which could be anything as it is not used by the code

afnTM: Transformation matrix in the format as given to affine2d [x1,y1,0;x2,y2,0;x3,y3,1] where x is the column number and y is row number( see https://in.mathworks.com/help/images/image-coordinate-systems.html for details)

output: name of the output image

toshow :  0 if don't want to show the output images generated
	  any other value then output images will be displayed in subplot


Syntax:(Examples)
For angle:
	A2a_Paras_2018CSB1111_2020_CS517('mozart.jpg',1,30,[],[],'test',1)
	A2a_Paras_2018CSB1111_2020_CS517('test1.jfif',1,30,[],[],'test',1)

For any affine matrix
	A2a_Paras_2018CSB1111_2020_CS517('mozart.jpg',2,[],[],[2,1,0;1,2,0;0,0,1],'test',1)
	A2a_Paras_2018CSB1111_2020_CS517('test1.jfif',2,[],[],[2,1,0;1,2,0;0,0,1],'test',1)
or
	A2a_Paras_2018CSB1111_2020_CS517('mozart.jpg',2,[],[],mat,'test',1)
	where mat is the affine matrix

Description of Method:
	I have used the affine2d Method for the transformation using the matrix [lambda*(c1-1)+1,lambda*c2,0;lambda*(c3),lambda*(c4-1)+1,0;lambda*c5,lambda*c6,1] varying the parameter from 0 to 1 in steps to get 
	different numbers of intermediate images. The number of images can be specified in the code by as the top which is currently set to 40, so lambda takes steps of 1/40=0.025 .
	lambda =0 the matrix becomes eye(3) which when applied gives input image only. when lambda becomes 1 it becomes [c1 ,c2 ,0;c3,c4,0;c5,c6,1]
	which is just the completely transformed image.Padding is providing in the image so that gif gets a stable background.
	For creating the gif I have used imwrite function of matlab and provided the 10 seconds per frame. It could be changed by easily modifying the parameter in the function.
	For better effects in gif I have append the images first in forward order and then in backward order which could be easily modified by comment a small paragraph in the code.


	