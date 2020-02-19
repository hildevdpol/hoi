%% *HC 3.3 & 3.4: TRANSFORMATIONS AND INTERPOLATION (IMAGE PROCESSING COURSE LAB EXERCISES)*


close all;
clear all;
clc;

%% Image Interpolation
% create a synthetic image with MATLAB's function |peaks|
[X,Y] = meshgrid(-3:3);
V = peaks(X,Y);
figure
surf(X,Y,V)
title('Original Sampling');

%%
% create a finer grid of 25x25
% pixels (in the interval [-3,3]x[-3,3]) and estimate the function/image 
% values for those locations.
[Xq,Yq] = meshgrid(-3:0.25:3);
Vq = interp2(X,Y,V,Xq,Yq);
figure
surf(Xq,Yq,Vq);
title('Cubic Interpolation Over Finer Grid');

%%
% The following illustration shows the placement of interpolated values 
% (in red) among nine original sample values (in black)
imshow(imread('trans_ntimes.png'));

%%
% Now, let us load a real image and study different interpolation methods. 
I = im2double(imread('trans_testimage.png'));
imshow(I,[]);
title('Test Image');

%%
%  return the interpolated values (image B) on a refined grid formed by repeatedly halving the 
% intervals k times in each dimension.
I_interpolated = interp2(I, 3);
imshow(I_interpolated,[]);
title('Linear Interpolation');



% create a simple 5x5 pixel test image, which only contains one
% bright pixel in the middle
I_test = zeros(5,5);
I_test(3,3) = 1;
imshow(I_test, 'InitialMagnification', 2000);

%%
% Now, apply the interpolation by inserting 7 interpolated points between the original sample values
figure;
I_linear = interp2(I_test, 3, 'linear');
subplot(2,2,1); imshow(I_linear,[]); title('Linear');

I_nearest = interp2(I_test, 3, 'nearest');
subplot(2,2,2); imshow(I_nearest,[]); title('Nearest');

I_cubic = interp2(I_test, 3, 'cubic');
subplot(2,2,3); imshow(I_cubic,[]); title('Cubic');

I_spline = interp2(I_test, 3, 'spline');
subplot(2,2,4); imshow(I_spline,[]); title('Spline');

%%
% This can also be done using the meshgrids, for example:
[X,Y] = meshgrid(1:5);
[Xq,Yq] = meshgrid(1:0.125:5);

I_linear  = interp2(X, Y, I_test, Xq, Yq, 'linear');
I_nearest = interp2(X, Y, I_test, Xq, Yq, 'nearest');

figure;
subplot(2,2,1); imshow(I_linear, []); title('Linear');
subplot(2,2,2); imshow(I_nearest, []); title('Nearest');
subplot(2,2,3); surf(Xq, Yq, I_linear); title('Linear');
subplot(2,2,4); surf(Xq, Yq, I_nearest); title('Nearest');


%% Image Transformations
% The basic image transformations include image translation, scaling and
% rotation. 
I_translated_nearest = imtranslate(I, [10,30.5], 'nearest');
I_translated_cubic =  imtranslate(I, [10,30.5], 'bicubic');
figure;
subplot(1,2,1); imshow(I_translated_nearest,[]); title('Translated, Nearest');
subplot(1,2,2); imshow(I_translated_cubic,[]); title('Translated, Cubic');


%%
% To scale the image with a certain factor
scale = 3;
I_scaled_nearest = imresize(I, scale, 'nearest');
I_scaled_cubic = imresize(I, scale, 'bicubic');
figure;
subplot(1,2,1); imshow(I_scaled_nearest,[]); title('Scaled 300%, Nearest');
subplot(1,2,2); imshow(I_scaled_cubic,[]); title('Scaled 300%, Cubic');

%%
% For image rotation: 
angle = 30;
I_rotated_nearest = imrotate(I, angle, 'nearest');
I_rotated_cubic =  imrotate(I, angle, 'bicubic');
figure;
subplot(1,2,1); imshow(I_rotated_nearest,[]); title('Rotated 30 degrees, Nearest');
subplot(1,2,2); imshow(I_rotated_cubic,[]); title('Rotated 30 degrees, Bicubic');

%%
% perform the scaling and rotation operations using |imwarp|:
scale = 2;
A_scale = [scale 0 0; 0 scale 0; 0 0 1]
tform = affine2d(A_scale);
I_tform_scale = imwarp(I, tform);
figure; imshow(I_tform_scale); title('Scaled image');

angle = 30 * pi / 180;
A_rotate = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1]
tform = affine2d(A_rotate);
I_tform_rotate = imwarp(I, tform);
figure; imshow(I_tform_rotate); title('Rotated 30 degrees image');

%%
% We can even combine those two operations in one:
A_combined = A_rotate * A_scale
tform = affine2d(A_combined);
I_tform = imwarp(I, tform);
figure; imshow(I_tform); title('Scaled and rotated image');

%%
% _*Exercise 3*: Interpolating a MRI scan._
%
% * _Find an MRI image of the liver (e.g using Google)._
% * _Apply 4 types of interpolation to produce 4 images of size of the original image 

liver1 = im2double(imread('lever.jpeg'));
liver=liver1(1:352, 1:500);
liver = padarray(liver, [74,0]);
figure, imshow(liver), title('MRI lever van Google');
liver_50 = imresize(liver, [50 50], 'nearest');
figure, imshow(liver_50), title('Lever 50x50');

step = 1/(510/50);
[X,Y] = meshgrid(1:1:50);
[Xq,Yq] = meshgrid(1:step:50);

liver_linear  = interp2(X, Y, liver_50, Xq, Yq, 'linear');
liver_nearest = interp2(X, Y, liver_50, Xq, Yq, 'nearest');
liver_cubic = interp2(X, Y, liver_50, Xq, Yq, 'cubic');
liver_spline = interp2(X, Y, liver_50, Xq, Yq, 'spline');

figure;
subplot(2,2,1); imshow(liver_linear, []); title('linear');
subplot(2,2,2); imshow(liver_nearest, []); title('nearest');
subplot(2,2,3); imshow(liver_cubic, []); title('cubic');
subplot(2,2,4); imshow(liver_spline, []); title('spline');

liver_diff1 = liver - liver_linear;
liver_diff2 = liver - liver_nearest;
liver_diff3 = liver - liver_cubic;
liver_diff4 = liver - liver_spline;

figure;
subplot(2,2,1); imshow(liver_diff1, []); title('Difference original and linear');
subplot(2,2,2); imshow(liver_diff2, []); title('Difference original and nearest');
subplot(2,2,3); imshow(liver_diff3, []); title('Differnce original and cubic');
subplot(2,2,4); imshow(liver_diff4, []); title('Differnce original and spline');

%%
% _*Exercise 4*: Rotating a MRI scan._
clear all; close all; clc; 
MRI1 = im2double(imread('lever.jpeg'));
MRI=MRI1(1:352, 1:500);
MRI = padarray(MRI, [74,0]);


% * _Smooth the image using the Gaussian blur filter imgaussfilt(A, 4)._
MRI_gausfilt = imgaussfilt(MRI, 4);
% * _Find the location of the pixel with the highest intensity ( choose any one of them._
[X,Y]=find(MRI_gausfilt==max(MRI_gausfilt(:)));
[maxintens,idmax] = max(MRI_gausfilt(:));
[rowmax, colmax] = find(MRI_gausfilt == maxintens);
% * _Translate the image using_ 
% * _Now rotate the image using counterclockwise for 45
% degrees and display the result._
value_max=MRI_gausfilt(322,15);

xtr = (500/2)-X;
ytr = (500/2)-Y;
translation = [xtr, ytr];
transMRI = imtranslate(MRI_gausfilt,translation);
I_liver_translated = imtranslate(MRI_gausfilt, [xtr,ytr], 'nearest');
 

I_liver_rotated=imrotate(I_liver_translated,45);

subplot(1,3,1); imshow(MRI_gausfilt), title('Gaussian blur filter');
subplot(1,3,2); imshow(I_liver_translated); title('Translatie met max waarde (322,15)');
subplot (1,3,3); imshow(I_liver_rotated); title('rotation');



%%
% _*Exercise 5*: Warping a MRI scan._
%

clear all
I1 = im2double(imread('lever.jpeg'));

I=I1(1:352, 1:500);
I = padarray(I, [74,0]);
% * _Scale the image to the size of 100x100 pixels and call it_ 
% 
I_orig = imresize(I, [100 100], 'nearest');
figure; imshow(I_orig);
% * _Rotate the original image
I_90 = imrotate(I_orig, 90, 'bilinear');


I1 = imrotate(I_orig,15,'bilinear');
I2 = imrotate(I1,15,'bilinear');
I3 = imrotate(I2,15,'bilinear');
I4 = imrotate(I3,15,'bilinear');
I5 = imrotate(I4,15,'bilinear');
I_90_6steps = imrotate(I5,15,'bilinear');
I_90_6steps = I_90_6steps(122:221, 122:221);
 

diff = I_90 - I_90_6steps;
figure;

subplot(1,3,1); imshow(I_90,'InitialMagnification', 2000), title('I 90');
subplot(1,3,2); imshow(I_90_6steps,'InitialMagnification', 2000), title('I 90 6steps');
subplot(1,3,3); imshow(diff,'InitialMagnification', 2000), title('difference');


%% * _Do the same procedure but with different interpolation scheme,_
clear all; close all; clc; 
I1 = im2double(imread('lever.jpeg'));

I=I1(1:352, 1:500);
I = padarray(I, [74,0]);
% * _Scale the image to the size of 100x100 pixels
I_orig = imresize(I, [100 100], 'nearest');
% * _Rotate the original image
I_90 = imrotate(I_orig, 90, 'bicubic');

I1 = imrotate(I_orig,15,'bicubic');
I2 = imrotate(I1,15,'bicubic');
I3 = imrotate(I2,15,'bicubic');
I4 = imrotate(I3,15,'bicubic');
I5 = imrotate(I4,15,'bicubic');
I_90_6steps = imrotate(I5,15,'bicubic');
I_90_6steps = I_90_6steps(122:221, 122:221);

diff = I_90 - I_90_6steps;

figure;
subplot(1,3,1); imshow(I_90,'InitialMagnification', 2000), title('I 90 bicubic');
subplot(1,3,2), imshow(I_90_6steps,'InitialMagnification', 2000), title('I 90 6steps bicubic');
subplot(1,3,3), imshow(diff,'InitialMagnification', 2000), title('Difference bicubic');