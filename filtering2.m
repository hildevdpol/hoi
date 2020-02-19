%% *2 FILTERING (IMAGE PROCESSING COURSE LAB EXERCISES)*
% assignment Hilde van der Pol, week 2

% 2.1 Loading images
% Clean up your environment:
clear all   % clears all variables from the workspace
close all   % closes all figures
clc         % clears the command window (remove/supress this command if you want to see the output in the command window of previous runs of the script)

%
% Read an example 2D image, which is provided by the IPT
I = imread('testpat1.png');

% Let's print a part of I to the workspace to check it's values.
I(10:15, 10:15)

%  convert |I| to a variable of class |double|
I = im2double(I);
I(10:15, 10:15)

% To visualize the image 
imshow(I);
title('Test image');

% 2.2 Gaussian filter
% apply a Gaussian blur filter to the loaded image use the function 
sigma = 2;
I_filter_Gauss = imgaussfilt(I, sigma);
figure;
imshow(I_filter_Gauss);
title('Filtered Image (Gaussian) sigma=2');


% 2.3 Median filter
I_filter_median = medfilt2(I, [3, 3]);
figure;
imshow(I_filter_median);
title('Filtered Image (Median)'); 


% 2.3 Comparison of Gaussian and Median filters
bin_bounds= 0:1; 
figure; 
subplot(1,2,1), histogram(I_filter_Gauss)
subplot(1,2,2), histogram(I_filter_median)
%gaussion filter blurred meer in het midden dan de median filter 


% _*Exercise 4:* Subtract the original image from both filtered images.
subG= I_filter_Gauss-I; 
subm= I_filter_median-I;
% subplot(1,2,1), imshow(subG)
% subplot(1,2,2), imshow(subm)


% _*Exercise 5:* Create two 'noisy' copies of image I; one by adding Gaussian noise
% with mean=0 and variance=0.1
% and the second by adding "salt and pepper" noise.  

InoiseG = imnoise(I,'gaussian',0,0.1);
Inoisem = imnoise(I,'salt & pepper',0.1);
% subplot(1,2,1), imshow(InoiseG)
% subplot(1,2,2), imshow(Inoisem)

%
% _*Exercise 6:* For both noisy images, apply a Gaussian and
% median filter (separately).
InoiseGgauss= imgaussfilt(InoiseG,sigma);
InoiseGmedian= medfilt2(InoiseG, [3, 3]);
Inoisemgauss= imgaussfilt(Inoisem,sigma);
Inoisemmedian= medfilt2(Inoisem, [3, 3]);

% figure; 
subplot(1,4,1), imshow(InoiseGgauss);
subplot(1,4,2), imshow(InoiseGmedian);
subplot(1,4,3), imshow(Inoisemgauss);
subplot(1,4,4), imshow(Inoisemmedian);

figure;
subplot(1,4,1), histogram(InoiseGgauss);
subplot(1,4,2), histogram(InoiseGmedian);
subplot(1,4,3), histogram(Inoisemgauss);
subplot(1,4,4), histogram(Inoisemmedian);

% 2.4 Filtering using user-defined kernels
h = fspecial('average', [5 5])

% In order to use this kernel for convolution with an image, 
% the function |B = imfilter(A,h)| can be used:
I_filter_h = imfilter(I, h);
figure;
imshow(I_filter_h);
title('Filtered Image (Custom average kernel)');

%
% In the same manner, we can get the kernel of the Gaussian filter that was used before:
h_Gauss = fspecial('gaussian', [5 5])
I_filter_h_Gauss = imfilter(I, h_Gauss);
figure;
imshow(I_filter_h_Gauss);
title('Filtered Image (Custom Gaussian kernel)');

% To visualize the kernel in 3D (as a surface) we can run:
figure;
surf(h_Gauss);
title('Gaussian kernel');

%
% We can also define our own kernels, for example to compute the image derivative in x-direction we 
% define kernel |h| as
h_dx = [-1 0 1; -1 0 1; -1 0 1]

%
% and to find the y-derivative, the kernel |h| is given as
h_dy = [-1 -1 -1; 0 0 0; 1 1 1]

%
% Let us calculate the image derivatives and the gradient 
I_dx = imfilter(I, h_dx);
I_dy = imfilter(I, h_dy);
I_gradient = sqrt(I_dx.^2+I_dy.^2);

%
% _*Question:* What is the largest value in the_ |I_dx| _or_ |I_dy| _that one can get assuming that the input image is double with values between 0 and 1?
maxI_dx= max(max(I_dx));
maxI_dy= max (max(I_dy));
%

figure;
subplot(1,3,1); imshow(I_dx,[]); title('x-derivative');
subplot(1,3,2); imshow(I_dy,[]); title('y-derivative');
subplot(1,3,3); imshow(I_gradient,[]); title('gradient magnitude');

%% 
% Similarly, we can also define a Gaussian kernel ourselves
sigma = 0.5;
sz = 2; % integer which defines the range of x and y in the form (-sz, sz)
[x,y] = meshgrid(-sz:sz,-sz:sz)  % get all possible combinations of (x,y) coordinates of the 2D kernel
h_gauss = (1/(2*pi*(sigma^2))) * exp(-(x.^2+y.^2)/(2*sigma^2))

% Create a surface plot of the Gaussian kernel
figure; surf(h_gauss); title('Gaussian kernel')

%
% _*Exercise 7:* Likewise, we can define the first Gaussian derivatives in both the x and y
% directions using the analytical expression for the Gaussian kernel:_
h_dx_gauss = -(x/(2*pi*sigma^4)) .* exp(-(x.^2+y.^2)/(2*sigma^2));
h_dy_gauss = -(y/(2*pi*sigma^4)) .* exp(-(x.^2+y.^2)/(2*sigma^2));

% - _Visualize the derivative kernels in 3D plots, use the same meshgrid as
% above._
subplot(1,2,1), surf(h_dx_gauss)
subplot(1,2,2), surf(h_dy_gauss)
%
% - _Convolve the original image with two derivative kernels and compute the
% gradient image._
I_gauss_x = imfilter(I, h_dx_gauss);
I_gauss_y = imfilter(I, h_dy_gauss);
I_gradientgauss = sqrt(I_gauss_x.^2+I_gauss_y.^2);

%
% - _Display the result and compare with the gradient image obtained using_
% |h_dx| _and_ |h_dy| _(see above)._
figure;
imshow(I_gradientgauss, [])

%% 2.5 Padding and Boundary Conditions

% Let us take a small crop of the original image and pad it using all
% mentioned ways
Ii = I(30:200,30:200);
Ii_zero = padarray(Ii, [15, 30]);
Ii_circ = padarray(Ii, [15, 30],'circular');
Ii_rep = padarray(Ii, [15, 30],'replicate');
Ii_sym = padarray(Ii, [15, 30],'symmetric');

figure;
subplot(3,2,1); imshow(Ii); title('Cropped image');
subplot(3,2,2); imshow(Ii_zero); title('Padding: zero');
subplot(3,2,3); imshow(Ii_circ); title('Padding: circular');
subplot(3,2,4); imshow(Ii_rep); title('Padding: replicate');
subplot(3,2,5); imshow(Ii_sym); title('Padding: symmetric');
%

%
% - _Crop the original image_ |I| _to_ [30, 200]x[30, 200] (_call it, for

Ii = I(30:200,30:200);
Ii_res = imfilter(Ii, h_dx, 'circular');
figure;
imshow(Ii_res, []); title('x-derivative with circular boundary conditions');
%
% _*Exercise 8*: Study the artifacts caused by different padding schemes. 
% For each of four padding types do the following:_

Ii = I(30:200,30:200);
% gaus filter 15x15 sigma 10
sigma = 10;
sz = 15; % integer which defines the range of x and y in the form (-sz, sz)
[x,y] = meshgrid(-sz:sz,-sz:sz);  % get all possible combinations of (x,y) coordinates of the 2D kernel
h_gauss = (1/(2*pi*(sigma^2))) * exp(-(x.^2+y.^2)/(2*sigma^2));

Ii_filtered = imfilter(Ii,h_gauss);

figure, surf(h_gauss);
figure, imshow(Ii_filtered), title('gaus filter 15x15 sigma 10');
%%
% - _Filter_ |Ii| _using the Gaussain filter with_ |sigma=10| _and the size of the filter_ |hsize=15|. _This will result in image_ |Ii_filtered| 
I_gausfilt = imfilter(I,h_gauss);
I_gold = I_gausfilt(30:200,30:200);
figure;
subplot(1,3,1), imshow(I_gold,[]), title('I gold');
I_diff = Ii_filtered - I_gold;
subplot(1,3,2), imshow(I_diff,[]), title('I diff');
subplot(1,3,3), imshow(Ii_filtered,[]), title('Ii');


%%
% - _For each of four padding types, present 2 images: 1) the result of
% filtering and 2) the difference image_ I_diff _with artifacts. In total, 8 plots._
Ii = I(30:200,30:200);
Ii_uitknipmetzeropadding = imfilter(Ii,h_gauss);
Ii_uitknipmetcircpadding = imfilter(Ii,h_gauss,'circular');
Ii_uitknipmetreppadding = imfilter(Ii,h_gauss,'replicate');
Ii_uitknipmetsympadding = imfilter(Ii,h_gauss,'symmetric');

figure;
subplot(2,4,1); imshow(Ii_uitknipmetzeropadding,[]); title('uitgeknipt en dan filter Padding: zero filtered');
subplot(2,4,2); imshow(Ii_uitknipmetcircpadding,[]); title('uitgeknipt en dan filter Padding: circular filtered');
subplot(2,4,3); imshow(Ii_uitknipmetreppadding,[]); title('uitgeknipt en dan filter Padding: replicate filtered');
subplot(2,4,4); imshow(Ii_uitknipmetsympadding,[]); title('uitgeknipt en dan filter Padding: symmetric filtered');

Ii_diffzero = Ii_uitknipmetzeropadding - I_gold;
Ii_diffcirc = Ii_uitknipmetcircpadding - I_gold;
Ii_diffrep = Ii_uitknipmetreppadding - I_gold;
Ii_diffsym = Ii_uitknipmetsympadding - I_gold;

subplot(2,4,5); imshow(Ii_diffzero,[]); title('Padding: zero difference');
subplot(2,4,6); imshow(Ii_diffcirc,[]); title('Padding: circular difference');
subplot(2,4,7); imshow(Ii_diffrep,[]); title('Padding: replicate difference');
subplot(2,4,8); imshow(Ii_diffsym,[]); title('Padding: symmetric difference');



%% 2.6 Real Applications (MRI with RF spikes)
% Now, let us study how the filtering can be used to estimate the
% radiofrequency (RF) pattern in a real MRI image. 
image_mri = imread('rfspike.jpg');
image_mri = im2double(image_mri(:,:,1));
figure; imshow(image_mri); title('RF spike MRI');

%
% _*Exercise 9:* Estimate the stripe-like noisy pattern from the image using the
% median filter by applying an appropriate_ |h|. 

figure;imshow(imread('rfspike_bg.bmp')); title('MRI estimated pattern');
figure;imshow(imread('rfspike_clean.bmp')); title('MRI with removed pattern');

EstpatternMRI = medfilt2(image_mri, [100, 1]);
MRIgoed = image_mri - EstpatternMRI;
figure;imshow(adapthisteq(EstpatternMRI)), title('MRI patroon');
figure;imshow(adapthisteq(MRIgoed)); title('MRI goed');

%% 2.7 Real Applications (MRI image noise/filtering)
% _*Exercise 10:* Find an MRI image of the liver (e.g using Google). 
picliver= imread('levergoed.jpeg');
liver=im2double(picliver); 
figure; 
imshow(liver), title 'Lever van google'

%
% _*Exercise 11:* Create three copies of the image with the following types of
% noise: 1) 'salt & pepper', 2) 'gaussian' and 3) 'motion'. 

liversalt = imnoise(liver,'salt & pepper', 0.05);
livergaus = imnoise(liver,'gaussian',0,0.01);
h = fspecial('motion',9,0);
livermotion = imfilter(liver,h);

figure;
subplot(1,3,1), imshow(liversalt,[]), title('salt&pepper noise');
subplot(1,3,2), imshow(livergaus,[]), title('gaussian noise');
subplot(1,3,3), imshow(livermotion,[]), title('motion noise');


%
% _*Exercise 12:*
% Apply median, Gaussian, Gaussian-x-derivative and Gaussian-y-derivative 
% filtering, to each of the three distorted images. 
a = 5;
b = a;

liversaltmed = medfilt2(liversalt,[a, b]);
livergausmed = medfilt2(livergaus,[a, b]);
livermotionmed = medfilt2(livermotion,[a, b]);

figure;
subplot(4,3,1), imshow(liversaltmed,[]), title('salt&pepper met median filter');
subplot(4,3,2), imshow(livergausmed,[]), title('gauss met median filter');
subplot(4,3,3), imshow(livermotionmed,[]), title('motion met median filter');

size = 3;
sigma = 1;
gausfilter = fspecial('gaussian',size,sigma);
liversaltgaus = imfilter(liversalt,gausfilter);
livergausgaus = imfilter(livergaus,gausfilter);
livermotiongaus = imfilter(livermotion,gausfilter);

subplot(4,3,4), imshow(liversaltgaus,[]), title('salt&pepper met gaussian filter');
subplot(4,3,5), imshow(livergausgaus,[]), title('gauss met gaussian filter');
subplot(4,3,6), imshow(livermotiongaus,[]), title('motion met gaussian filter');

h_dx_gaus = -(x/(2*pi*sigma^4)) .* exp(-(x.^2+y.^2)/(2*sigma^2));

liversaltgausxder = imfilter(liversalt,h_dx_gaus);
livergausgausxder = imfilter(livergaus,h_dx_gaus);
livermotiongausxder = imfilter(livermotion,h_dx_gaus);

subplot(4,3,7), imshow(liversaltgausxder,[]), title('salt&pepper met gaussian X der filter');
subplot(4,3,8), imshow(livergausgausxder,[]), title('gauss met gaussian X der filter');
subplot(4,3,9), imshow(livermotiongausxder,[]), title('motion met gaussian X der filter');

h_dy_gaus = -(y/(2*pi*sigma^4)) .* exp(-(x.^2+y.^2)/(2*sigma^2));

liversaltgausyder = imfilter(liversalt,h_dy_gaus);
livergausgausyder = imfilter(livergaus,h_dy_gaus);
livermotiongausyder = imfilter(livermotion,h_dy_gaus);

subplot(4,3,10), imshow(liversaltgausyder,[]), title('salt&pepper met gaussian Y der filter');
subplot(4,3,11), imshow(livergausgausyder,[]), title('gauss met gaussian Y der filter');
subplot(4,3,12), imshow(livermotiongausyder,[]), title('motion met gaussian Y der filter');
