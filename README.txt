This software package titled: JISR, performs 2D and 3D joint image segmentation and registration for a pair of images, either synthetic or medical, which are of the same size. The user should put the images in the images 2D and images3D folder while running the code. 

The code is tested with the MATLAB 2018a version. Some functions use parallelization to accelerate computation. This is done by using "parfor" loops. Also, we have use mex functions that are automatically generated from the .m file using the Matlab Coder toolbox. The compiled .mex files and the original .m files are both provided. 

General instructions to run the code:
1. M = Moving image and F = Fixed image, phi = Level set function defining the contour of the segmeted source image
2. The default number of maximum refinement levels are set as 3. 
3. To run the files, run the appropriate mainfunction_*.m files for each of the example listed below.
4. For the 3D images: If the compiled .mex files do not work with your Matlab version or achitecture, the compilemex.m file in the run/ folder can be run. This function compiles all the MEX functions using the OpenMP parallel framework, provided that a supported compiler is installed. If compilation fails, the native Matlab .m files can be used by removing the "_mex" suffix from the function names.

For different examples, we have used different main files for running with the respective settings:

1. Circle to White C example (Fig. 2 in article)
run file: mainfunction_whitec.m function is used to run this example. The image files are in images2D folder. The intensities are normalized. The error tolerance value is set as tol = 1e-04 and the maximum iterations are set as 50. Users can set it accordingly on lines 342 onwards in the MultiResolution2D_whitec.m. The rest of the parameters are set in the file setparamters_whitec.m script in the setparameters folder. 


2. Circle to Triangle example (Fig. 3 & 4 in article):
run file: mainfunction_triangle.m function is used to run this example. The image files are in images2D folder. The images with tdifferent noise intensities are given in .mat files. The images are given as triangle_*.mat files. Modify the source and target images in the mainfunction to the correct file names. The intensities are normalized. The error tolerance value is set as tol = 1e-04 and the maximum iterations are set as 50. Users can set it accordingly on lines 342 onwards in the MultiResolution2D_whitec.m. To set the parameters, modify the setparameters_triangle.m function. 

3. Brain MRI with multiple sclerosis taken from Brainweb website (Fig. 5 in article)
run file: mainfunction_brainweb.m function is used to run this example. Download the files from the website http://brainweb.bic.mni.mcgill.ca/brainweb/anatomic_ms.html. Save the files in images2D folder. The discrete model is stored as _rawb. The white matter segmented image is also obtained for the corresponding subject. To set the parameters, modify the setparameters_brainweb.m function.

4. Lung CT image (Fig. 6 in article)
run file: mainfunction_lung.m function is used to run this example. The images are obtained at the inhale and exhale stage of the breathing cycle. Download the dataset of the lung images from the website: https://www.dir-lab.com/ReferenceData.html
Store all the files in images2D folder. The segmentation of the lung images for the evaluation of phi is done using level set segmentation method. The code for segmentation can be found in  http://www.imagecomputing.org/~cmli/. To set the parameters, modify the setparameters_lung.m function.


5. Sphere to Bunny Image (Fig. 7 in article)
run file: mainfunction_synthetic.m function is used to run this example. 
The image files are found in images3D folder. The output image matrices are then registered in the non-rigid framework in our code. To set the parameters, modify the setparameters_synthetic.m function.

6. Brainweb 3D brain MRI images (Fig. 7)
run file: mainfunction_medical.m function is used to run this example. The images are downloaded from http://brainweb.bic.mni.mcgill.ca/brainweb/anatomic_normal_20.html website. The input images and the gray matter segmented images are obtained from the website. Download the images in the _.rawb format. 
