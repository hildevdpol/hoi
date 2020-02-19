%% *3 BINAIRY MORPHOLOGY (IMAGE PROCESSING COURSE LAB EXERCISES)*


%% Create image
% Clean up your environment:
clear all
close all
clc

%%
% We first create a binary image with a square in the middle
I = false(250); % create a 250x250 logical matrix of zeros
I(50:200,50:200) = true;
I(100:150,100:150) = false;
figure; imshow(I)

%% Structuring elements
% Create a disk structuring element with |R=20| and |N=0| using:
se = strel('disk',20,0);

% Visualize the (neighborhood matrix of the) disk-shaped strel:
nhood_disk = getnhood(se);
figure; imshow(nhood_disk); title('disk-shaped strel')



%% Erosion and dilation
% Erode and dilatie image |I| using the disk-shaped strel with |R=20| and
% |N=0|
I_ero = imerode(I, se);
I_dil = imdilate(I, se);

%%
% Visualize the original, eroded and dilated images next to each other
figure;
subplot(1,3,1); imshow(I); title('original image')
subplot(1,3,2); imshow(I_ero); title('eroded image')
subplot(1,3,3); imshow(I_dil); title('dilated image')


%% Opening and closing
% Let us take another (grey-scaled) image and segment/binarize it 
% (threshold at 25% of the maximum image intensity)
Is = im2double(imread('snowflakes.png'))>0.25;
figure; imshow(Is)

%%
% Apply the 'opening' and 'closing' morphological operation to image |Is|
% using a disk-shaped structuring element with a radius of 2 pixels
se = strel('disk',2,0);
Is_open = imopen(Is,se);
Is_close = imclose(Is,se);

%%
% Visualize the results
figure;
subplot(3,1,1); imshow(Is); title('original image')
subplot(3,1,2); imshow(Is_open); title('opened image')
subplot(3,1,3); imshow(Is_close); title('closed image')


%% Boundary detection

% Load the image
Ic = imread('circles.png');

% Define a disk-shaped strel with a 1 pixel radius
se = strel('disk',1,0);

% Obtain an eroded and dilated version of image |Ic|
Ic_ero = imerode(Ic,se);
Ic_dil = imdilate(Ic,se);

% Obtain the thick and inner boundary
bound_thick = Ic_dil-Ic_ero;
bound_inner = Ic-Ic_ero;

% Visualize the boundaries
figure;
subplot(1,2,1); imshow(bound_thick); title('thick boundary')
subplot(1,2,2); imshow(bound_inner); title('inner boundary')




%% Real applications (liver cyst boundary)
% _*Exercise 8:*_
%
% * _Read the liver_cyst.png 
% * _Segment the cyst by thresholding 
% * _Fill holes in the images 
% * _Derive the outer boundary of the opened image 
% * _Create and visualize the original image with an overlay of
% boundary 
Iliver = imread('liver_cyst.png');
Iliver = im2double(Iliver);

Iliversegmented = logical(Iliver<0.15);
Iliverfill = imfill(Iliversegmented,'holes');

se = strel('disk',5,0); 
Iliveropen = imopen(Iliverfill,se);

Iliverdil = imdilate(Iliveropen,se);
Iboundary = (Iliverdil - Iliversegmented);
Iliverboundary = Iliver + Iboundary;

figure;
subplot(2,3,1); imshow(Iliver); title('original')
subplot(2,3,2); imshow(Iliversegmented); title('segmented ')
subplot(2,3,3); imshow(Iliverfill); title('filled')
subplot(2,3,4); imshow(Iliveropen); title('opened');
subplot(2,3,5); imshow(Iliverboundary); title('boundary & overlay image');

%% Real applications (angiogram skeleton)
% _*Exercise 9:* The following exercise aims to find the 'skeleton' of
% an X-ray coronary angiogram. 

sigmaGauss = 1;
binthr = 0.4;
radiusClose = 5;
radiusOpen = 2;

% Step 1: Load image of an angiogram and convert to double
Iangio = im2double(imread('angiogram.jpg'));

% Step 2: Crop the image
Iangio_crop = Iangio(540:774,109:315);

% Step 3: Smooth the image using a Gaussian filter
Iangio_gaus = imgaussfilt(Iangio_crop,sigmaGauss);

% Step 4: Segment the cyst by simple thresholding at 40% of the intensity range
Iangio_bin = Iangio_gaus<binthr;

% Step 5: Fill the image to remove gaps
Iangio_fill = imfill(Iangio_bin,'holes');

% Step 6: Close the image to connect different segments using a disk-shaped strel
se=strel('disk',radiusClose,0);
Iangio_close = imclose(Iangio_fill,se);

% Step 7: Open the image to remove spurious segments using a disk-shaped strel
se=strel('disk',radiusOpen,0);
Iangio_open = imopen(Iangio_close,se);

% Step 8: Thinning of the image until your left with a skeleton
Iangio_skel = bwmorph(Iangio_open,'thin',inf);

% Visualize
figure
subplot(2,4,1); imshow(Iangio); title('original')
subplot(2,4,2); imshow(Iangio_crop); title('cropped')
subplot(2,4,3); imshow(Iangio_gaus); title('filtered')
subplot(2,4,4); imshow(Iangio_bin); title('binarized')
subplot(2,4,5); imshow(Iangio_fill); title('filled')
subplot(2,4,6); imshow(Iangio_close); title('closed')
subplot(2,4,7); imshow(Iangio_open); title('opened')
subplot(2,4,8); imshow(Iangio_skel); title('skeleton')


%%
% _*Exercise 10*: Search for an image of another type of angiogram or venogram and find the
% skeleton of a small part of the image

sigmaGauss = 0.9;
binthr = 0.5;
radiusClose = 1;
radiusOpen = 1;

% Step 1: Load image of an angiogram and convert to double
IangioMRI = im2double(imread('leg.jpeg'));

% % Step 2: Crop the image
Iangio_crop = IangioMRI(50:125,400:500)

% Step 3: Smooth the image using a Gaussian filter
Iangio_gaus = imgaussfilt(Iangio_crop,sigmaGauss);

% Step 4: Segment the cyst by simple thresholding at 40% of the intensity range
Iangio_bin = Iangio_gaus>binthr;

% Step 5: Fill the image to remove gaps
Iangio_fill = imfill(Iangio_bin,'holes');

% Step 6: Close the image to connect different segments using a disk-shaped strel
se=strel('disk',radiusClose,0);
Iangio_close = imclose(Iangio_fill,se);

% Step 7: Open the image to remove spurious segments using a disk-shaped strel
se=strel('disk',radiusOpen,0);
Iangio_open = imopen(Iangio_close,se);

% Step 8: Thinning of the image until your left with a skeleton
Iangio_skel = bwmorph(Iangio_open,'thin',inf);

% Visualize
figure;
axis tight
subplot(2,4,1); imshow(IangioMRI); title('original')
subplot(2,4,2); imshow(Iangio_crop); title('cropped')
subplot(2,4,3); imshow(Iangio_gaus); title('smooth')
subplot(2,4,4); imshow(Iangio_bin); title('segmented')
subplot(2,4,5); imshow(Iangio_fill); title('filled')
subplot(2,4,6); imshow(Iangio_close); title('closed')
subplot(2,4,7); imshow(Iangio_open); title('opened')
subplot(2,4,8); imshow(Iangio_skel); title('thinning')
