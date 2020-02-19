%% *HC 3.1 & 3.2: SEGMENTATION (IMAGE PROCESSING COURSE LAB EXERCISES)*


% Clean your environment
clear all
close all
clc

%%
% Let us load an example image that contains objects of interest placed on
% non-uniform noisy background. 
I = imread('rice.png');
imshow(I)

% histogram of image intensities already shows us that there is no
% clear separation between two classes. 
figure;
histogram(I)
title('Histogram of Intensity Values');

%%
% In order to segment the objects, we will try to threshold the original image 
threshold_value = 0.5;
bw = im2bw(I, threshold_value);
figure; imshow(bw); title('Thresholded image');


%%
% The other option is to apply more advanced thresholding methods, which
% automatically estimate the best value for |threshold_value|
threshold_level = graythresh(I)
bw = im2bw(I, threshold_level);
figure; imshow(bw); title('Thresholded using Otsu');

%%
%still not perfect, therefor: 
% In order  to segment the objects properly, we first have to estimate the
% non-uniform background and subtract it from the original image. 
background = imopen(I,strel('disk', 15));
figure;
surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
ax = gca;
ax.YDir = 'reverse';

% Now, let us plot the original image, the estimated background and their
% difference:
figure;
subplot(1,3,1);imshow(I); title('Original image');
subplot(1,3,2);imshow(background); title('Estimated background');
subplot(1,3,3);imshow(I - background); title('Difference');

%% 
% In order to segment the objects of interest, we have to binarize the 
% difference image. 
bw = im2bw(I - background, 0.19);
figure; imshow(bw); title('Segmented objects');

%%
% The value |0.19| in the previous example could be found either manually or
% by using the Otsu's automatic thresholding method. 
level = graythresh(I - background)

%%
% Now, we can apply MATLAB's methods to count the objects  
cc = bwconncomp(bw, 4)

%%
% To compute the area occupied by each object:
graindata = regionprops(cc,'Area');
graindata(1).Area % area of object with id=1
graindata(10).Area % area of object with id=10

%%
% To explore the histogram with the Area values we can do the following:
figure;
grain_areas = [graindata.Area]; % convert to the standard MATLAB vector
histogram(grain_areas)
title('Histogram of Object Areas');

%%
% Next, let us visualize a specific object.  To display only that object:
grain = false(size(bw)); % creates empty binary image with the same size as bw 
id=24;
grain(cc.PixelIdxList{id}) = true; % mark all the pixels in grain that correspond to object with idx
figure; imshow(grain); title('Object with id=24');

%%
% _*Exercise 1*: Segmentation of cells in color images:_ 
%Extract the B (blue) color image (
clear all; close all; clc; 

I_cell= imread('segmentation_cells.jpg');
[rows_cell, cols_cell, nchn_cell] = size(I_cell); 
I_cell_blu= im2double(I_cell(1:rows_cell, 1:cols_cell, 3));


% * _Randomly select a sub-region
I_cell_crop= I_cell_blu(1:200, 1:200); 


% * _Find the threshold which separates the cells from the background as
% good as possible 
threshold_level = graythresh(I_cell_crop)
bw = im2bw(I_cell_crop, threshold_level);

% * _ there are still some 'touching' cells in the binary image that have to 
% be split in separate objects. Using morphological operations
SE= strel('disk', 4,0); 
I_res= imerode(bw,SE); 
figure; imshow(I_res); title('I res');

% * _Count the number of segmented objects (
cc = bwconncomp(bw, 4);
cc2 = bwconncomp(I_res, 4);
% * _Find the object with the maximum area and produce 
% a new image which contains only object(s) with that size_
 cell_data = regionprops(cc2,'Area');
 cell_areas = [cell_data.Area];
 [Y, I]=max(cell_areas);
 
I_max = false(size(bw)); % creates empty binary image with the same ...
    %size as bw
id=27;
I_max(cc.PixelIdxList{id}) = true; % mark all the pixels in grain ...
    %that correspond to object with idx
figure; imshow(I_max); title('Image max');
 
% plotting
figure;
subplot(2,2,1); imshow(I_cell_crop); title('Image cropped'); 
subplot(2,2,2); imshow(bw); title('bw'); 
subplot(2,2,3); imshow(I_res); title('Image res, number segm= 33'); 
subplot(2,2,4); imshow(I_max); title('Image max'); 

%%
% _*Exercise 2*: Brain structure segmentation_
%
% _determine the fraction of gray matter, white
% matter and CSF in a T1-weighted brain image._ 

brain=imread('brain.jpg');
brain=im2double(brain); 
figure; imshow(brain);
% * _Using thresholding, segment the image into gray matter, white matter 
% and CSF. 
I_brain=rgb2gray(brain);
figure; imshow(I_brain);
figure; histogram(I_brain);
CSF = I_brain < 0.41 & I_brain > 0.05;

grey = I_brain < 0.58 & I_brain > 0.43;
white = I_brain < 0.8 & I_brain > 0.6;

% * _Determine their relative fractions with respect to the total
% brain volume 
fr_CSF=sum(CSF)/(sum(CSF)+sum(white)+sum(grey));
fr_White=sum(white)/(sum(CSF)+sum(white)+sum(grey));
fr_Grey=sum(grey)/(sum(CSF)+sum(white)+sum(grey));

fr_tot=fr_CSF+fr_White+fr_Grey;

% * _Copy the original and the three segmented images to Powerpoint. Put the
% relative fractions of these segments in the figure titles._
figure(3);
subplot(2,2,1); imshow(I_brain); title('original');
subplot(2,2,2); imshow(CSF); title('segmented CSF, fraction =0.15');
subplot(2,2,3); imshow(grey); title('segmented grey matter, fraction= 0.37');
subplot(2,2,4); imshow(white); title('segmented white matter, fraction=0.48');

%%
% _*Exercise 3*: Lymphnode segmentation in a histological image._
%
% _The final purpose of this exercise is to segment lymphnode tissue in a 
% histological image and obtain characteristics of this tissue.
clear all; close all; clc; 

he= imread('he_lympnode.jpg'); 

% * _Crop the image to remove the scaling-bar. convert the image to 
% first a grayscale image, than a |double| image._
he_lymp= he(1:1000, 1:1163); 
he_bw=im2bw(he_lymp); 
he_lympn=im2double(he_lymp);



% * _Segment all nuclei 
I_nuclei=he_lympn < 0.68; 

% * _Segment the lympnodes, i.e. create |I_lymphnode|. 
I_lympnode = he_lympn < 0.98 & he_lympn > 0.8;
threshold_level = graythresh(he_lympn);
I_lympnodeo = im2bw(he_lympn, threshold_level);

%erode to lose dots
se = strel('disk',4,0);
I_lympnode_ero = imerode(I_lympnodeo, se);

%preparing for filling in the next step  
se = strel('disk',5,0);
I_lympnode_closing1=imclose(I_lympnode_ero,se);

%fill in holes  
I_lympnode_fill=imfill(I_lympnode_closing1, 'holes');

%lose the little dots  
se = strel('disk',6,0);
I_lympnode_open=imopen(I_lympnode_fill,se);

%filling in the lymphnode  
se = strel('disk',13,0);
I_lympnode_closing=imclose(I_lympnode_open,se);

figure(5); 
subplot(3,3,1); imshow(I_lympnode); title('lympnode thresholded manual')
subplot(3,3,2); imshow(I_lympnodeo); title('lympnode thresholded Otsu')
subplot(3,3,3); imshow(I_lympnode_ero); title('lympnode using erosie');
subplot(3,3,4); imshow(I_lympnode_closing1); title('lympnode using closing');
subplot(3,3,5); imshow(I_lympnode_fill); title('lympnode using filling'); 
subplot(3,3,6); imshow(I_lympnode_open); title('lympnode using opening'); 
subplot(3,3,7); imshow(I_lympnode_closing); title('lympnode using closing');

% * determine the total area of% the lymphnode tissue and the areas of the two largest lymphnodes. 

%area: 
cc = bwconncomp(I_lympnode_closing, 4);
% * _Find the object with the maximum area and produce 
% a new image_ |I_max| _which contains only object(s) with that size_
graindata = regionprops(cc,'Area');
grainareas = [graindata.Area];

[gesorteerde_areas, idnum] = sort(grainareas);


% * _Use |I_lympnode| and  |I_nuclei| to create |I_nuclei_tissue| and |I_nuclei_lymphnode|. 

drielymfen = idnum(11:13);
grain = false(size(I_lympnode_closing));

grain1 = false(size(he_lympn));  
grain2 = false(size(he_lympn)); 
grain3 = false(size(he_lympn)); 

grain1(cc.PixelIdxList{drielymfen(1)}) = true;
grain2(cc.PixelIdxList{drielymfen(2)}) = true;
grain3(cc.PixelIdxList{drielymfen(3)}) = true;


figure; imshow(grain1+grain2+grain3); title('3 lymfeknopen');


I_nuclei_tissue = abs(he_lympn - I_lympnode_closing);
I_tissue=I_nuclei_tissue - I_nuclei;
I_nuclei_lympnode= abs(he_lympn-(~I_nuclei_tissue)-I_tissue);
figure; 
subplot(1,2,1); imshow(I_nuclei_tissue); title('image of nucleii with tissue');
subplot(1,2,2); imshow(I_nuclei_lympnode); title('image of nucleii with lympnode');

figure; histogram(I_nuclei_tissue); title('histogram of image with nucleii and tissue');

Cyto1= (0.75 < I_nuclei_tissue);
Cyto2=(I_nuclei_tissue<= 0.85); 
I_Cyto=Cyto1.*Cyto2;

Nuclei1= (0.55 < I_nuclei_lympnode);
Nuclei2=(I_nuclei_lympnode<= 0.65); 
I_Nuclei=Nuclei1.*Nuclei2;

figure;
subplot(1,2,1); imshow(I_Cyto); title('segmented nucleii normal tissue');
subplot(1,2,2); imshow(I_Nuclei); title('segmented nucleii lymphnode');


% * _Let us assume that a pixel in the original image either represents cytoplasm 
% or nuclei. Use the segmented images to calculate the nuclei-to-plasma ratio 
fr_cyto=sum(I_Cyto)/(sum(I_Cyto)+sum(I_Nuclei));
fr_nuclei=sum(I_Nuclei)/(sum(I_Cyto)+sum(I_Nuclei));


