%% exercise 1
clear all; close all; clc; 

% load the rectangle image
imgRect = imread('rectangle.png');
% display the image using imshow()
figure; imshow(imgRect);
 
% Fourier transform
specRect = fft2(imgRect);
% compute the magnitude of specRect
magRect = abs(specRect);

% color-scaled visualization
figure; imagesc(magRect);

% visualization with zero-frequency component in the center
figure; imagesc(fftshift(magRect));
figure; mesh(fftshift(magRect));

%% exercise 2 
% generate a Gaussian kernel
imgGauss = fspecial('gaussian', [5,5], 1);

% Fourier transform of the Gaussian kernel
specGauss = fft2(imgGauss, 512, 512);
magspecGauss= abs(specGauss);
figure; mesh(fftshift(magspecGauss));

%% Exercise 3 
% open the supplied interface
retainCoeffsGui


%% exercise 4

clear all; close all; clc; 
% create a spectrum with two peaks on the horizontal frequency axis
specH2 = zeros(512, 512);
specH2(257, 257 - 2) = 1;
specH2(257, 257 + 2) = 1;

% create a spectrum with two peaks on the vertical frequency axis
specV2 = zeros(512, 512);
specV2(257 - 2, 257) = 1;
specV2(257 + 2, 257) = 1

% compute inverse Fourier transform
imgH2 = ifft2(ifftshift(specH2));
imgV2 = ifft2(ifftshift(specV2));

% visualize the resulting images
figure;
subplot(2,2, 1); imagesc(imgH2); title('H2 image');
subplot(2,2, 2); mesh(imgH2); title('H2 mesh');
subplot(2,2, 3); imagesc(imgV2); title('V2 image');
subplot(2,2, 4); mesh(imgV2); title('V2 mesh');

%% ex 4 deel 2
specMoonMod = zeros(512, 512);
specMoonMod(257 + 2, 257 + 2) = 1;
specMoonMod(257 - 2, 257 - 2) = 1;

% compute inverse Fourier transform
imgDIA = ifft2(ifftshift(specMoonMod));

% visualize the resulting images
figure;
subplot(1,2, 1); imagesc(imgDIA); title('Diagonal imagesc');
subplot(1,2, 2); mesh(imgDIA); title('Diagonal mesh');

%% exercise 5
clear all; close all; clc; 

% Load brain.jpg and apply Gaussian noise
imgBrain = imread('brain.jpg');
imgBrain = imnoise(imgBrain, 'gaussian');

% Gaussian filtering in the spatial domain via convolution
filterGauss = fspecial('gaussian', [5,5], 1);
imgBrainConv = conv2(double(imgBrain), filterGauss, 'same');
imgBrainConv = uint8(imgBrainConv);


% Fourier transform of the noisy image
specBrain = fft2(imgBrain);
% Fourier transform of the Gaussian filter kernel
specFilter = fft2(filterGauss, 512, 512);

% Frequency-wise multiplication
specBrainMod = specBrain .* specFilter;

imgBrainMod = uint8(ifft2(specBrainMod));
figure; 
subplot (1,3,1); imshow(imgBrain);title ('Noisy')
subplot (1,3,2); imshow(imgBrainConv);title ('convolution')
subplot (1,3,3); imshow(imgBrainMod);title('Multiplication')

%% exercise 6
clear all; close all; clc; 
imgMoon = imread('moon.png');

% visualize log-scaled magnitude of Fourier spectrum
specMoon = fft2(imgMoon);
magMoon = abs(specMoon);
figure; imagesc(log10(fftshift(magMoon)));

% locatie noise peaks: (237, 237),(277, 237), (237, 277),  (277, 277) 
specMoonMod = fftshift(specMoon);
specMoonMod(237,237) = 4.2;
specMoonMod(237,277) = 4.2;
specMoonMod(277,237) = 4.2;
specMoonMod(277,277) = 4.2;

magMoonMod = abs(specMoonMod);
figure, imagesc(log10(fftshift(magMoonMod))), title('spectrum zonder ruis');

% inverse transform the modified spectrum
imgMoonMod = uint8(ifft2(ifftshift(specMoonMod)));
figure; imshow(imgMoonMod);

%% exercise 7
% write your code here
clear all; close all; clc; 

cat= imread('cat.jpg');
dog=imread('dog.jpg');

speccat= fft2(cat);
specdog=fft2(dog); 

% split spectrum into magnitude and phase
magcat = abs(speccat);
phasecat = angle(speccat);
magdog = abs(specdog);
phasedog = angle(specdog);

% reconstruct spectrum from magnitude and phase
spectrum1= magcat .* exp(i .* phasedog);
spectrum2= magdog .* exp(i .* phasecat);

img1 = uint8(real(ifft2(spectrum1)));
img2= uint8(real(ifft2(spectrum2)));

figure; 
subplot(2,2,1); imshow(cat); title('image cat')
subplot(2,2,2); imshow(dog); title('image dog')
subplot(2,2,3); imshow(img1); title('magnitude cat w/ phase dog')
subplot(2,2,4); imshow(img2); title('magnitude dog w/ phase cat')

