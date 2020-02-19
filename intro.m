%% hilde van der pol, assignment week 1

  % Reading a 2D image
  I_ileum = imread('images/practicum-histologie-3e-jaar-KT-rat_ileum_200x_HE.jpg');
  I_maag = imread('images/practicum-histologie-3e-jaar-KT-rat-maag-HE.jpg');
  I_prox= imread('images/practicum-histologie-3e-jaar-KT-rat-proximaal colon-HE.jpg');

  
   % dimensions of image
  num_dim_ileum = ndims(I_ileum);
  num_dim_maag = ndims(I_maag);
  num_dim_prox = ndims(I_prox);
 
  
% show image
figure(1);
subplot(1,3,1), imshow(I_ileum)
subplot(1,3,2), imshow(I_maag)
subplot(1,3,3), imshow(I_prox)

  % dimensions of image
  num_dim_ileum = ndims(I_ileum);
  num_dim_maag = ndims(I_maag);
  num_dim_prox = ndims(I_prox);
  
  % number of dimensions
disp(['number of dimensions ileum: ', num2str(num_dim_ileum)]);
% size of image
[rows_ileum, cols_ileum, nchn_ileum] = size(I_ileum);
disp(['rows: ', num2str(rows_ileum), ' columns: ', num2str(cols_ileum)]);
disp(['number of color channels: ', num2str(nchn_ileum)]);

  % number of dimensions
disp(['number of dimensions maag: ', num2str(num_dim_maag)]);
% size of image
[rows_maag, cols_maag, nchn_maag] = size(I_maag);
disp(['rows: ', num2str(rows_maag), ' columns: ', num2str(cols_maag)]);
disp(['number of color channels: ', num2str(nchn_maag)]);

  % number of dimensions
disp(['number of dimensions prox: ', num2str(num_dim_prox)]);
% size of image
[rows_prox, cols_prox, nchn_prox] = size(I_prox);
disp(['rows: ', num2str(rows_prox), ' columns: ', num2str(cols_prox)]);
disp(['number of color channels: ', num2str(nchn_prox)]);

% color channels
red_chn = I_ileum(1:rows_ileum, 1:cols_ileum, 1);
grn_chn = I_ileum(1:rows_ileum, 1:cols_ileum, 2);
blu_chn = I_ileum(1:rows_ileum, 1:cols_ileum, 3);
figure(2);
subplot(1, 3, 1), imshow(red_chn), title('Red Channel');
subplot(1, 3, 2), imshow(grn_chn), title('Green Channel');
subplot(1, 3, 3), imshow(blu_chn), title('Blue Channel');

% access the value of the matrix element at row: 471 and column: 599
value = red_chn(471, 599);
disp(value);

% take a region from rows 125 to 145 and colums 525:575
windw = red_chn(125:145, 525:575);
%also we can use: windw = red_chn(471-10:471+10, 599-10:599+10)
disp(windw);
figure, imshow(windw);

% %exercise 3 
I_ileum(size(I_ileum,1)/2-50:size(I_ileum,1)/2+50,size(I_ileum,2)/2-50:size(I_ileum,2)/2+50)=I_ileum(size(I_ileum,1)/2-50:size(I_ileum,1)/2+50,size(I_ileum,2)/2-50:size(I_ileum,2)/2+50)+100;imshow(I_ileum)

% %exercise 3 part 2  
 I_ileum(size(I_ileum,1)/2-50:size(I_ileum,1)/2+50,size(I_ileum,2)/2-50:size(I_ileum,2)/2+50)=I_ileum(size(I_ileum,1)/2-50:size(I_ileum,1)/2+50,size(I_ileum,2)/2-50:size(I_ileum,2)/2+50)-100;imshow(I_ileum)

red_thr = 170; % threshold for red channel
filt_im = uint8(zeros(size(I_ileum))); % initialize
filt_mask = uint8(red_chn >= red_thr); % creating our mask using our threshold
I_thr1 = I_ileum(:, :, 1).*filt_mask;  % element-wise multiplication on red chn
I_thr2 = I_ileum(:, :, 2).*filt_mask;  % element-wise multiplication on green chn
I_thr3 = I_ileum(:, :, 3).*filt_mask;  % element-wise multiplication on blue chn
filt_im(:,:,1) = I_thr1;
filt_im(:,:,2) = I_thr2;
filt_im(:,:,3) = I_thr3;
figure,
subplot(1,4,1), imshow(I_ileum); title('Ileum Image');
subplot(1,4,2), imshow(filt_im); title('Red Threshold');
grn_thr = 190;
filt_mask = uint8((red_chn>=red_thr).*(grn_chn<=grn_thr)); % thresholding green chn
I_thr1 = I_ileum(:, :, 1).*filt_mask;  % element-wise multiplication on red chn
I_thr2 = I_ileum(:, :, 2).*filt_mask;  % element-wise multiplication on green chn
I_thr3 = I_ileum(:, :, 3).*filt_mask;  % element-wise multiplication on blue chn
filt_im(:,:,1) = I_thr1;
filt_im(:,:,2) = I_thr2;
filt_im(:,:,3) = I_thr3;
subplot(1,4,3), imshow(filt_im); title('Red+Green Threshold');

% display filt_mask as a binary image
subplot(1,4,4); imshow(logical(filt_mask));
title('Binary Image of Red+Green Threshold');