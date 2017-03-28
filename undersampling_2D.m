clear all;
close all;
clc;

% Meghasyam Tummalacherla
% MRI Lab Project

% Undersampling and Compressed sensing in MRI

% In this code I will try to find the optimal 2D undersampling pattern
% Then reconstruct the brain from the undersampled data
% using the 2D constant undersampling pattern.

% INPUT - Reading all my brain image data
for i=1:172
    str=strcat('s6192/','*.MRDC.',num2str(i));
    listing(i) = dir(str);
end

fImage=zeros(256,256,length(listing)); % 3D data of my brain
for i=1:length(listing)
    fImage(:,:,i)=dicomread(strcat('s6192/',listing(i).name));
end

% To find the undersampling pattern U I will assume that all the images are
% of size 256x256, and have the same bandwidth so that when I calculate
% U,it can be used uniformly across the 3D data

figure(1);
% subplot(1,2,1);
imagesc(fImage(:,:,96));
title('Sample image to show the orientation of kx-ky space');
% subplot(1,2,2);
% imagesc(abs(ft(fImage(:,:,96))));
% xlabel('kx');
% ylabel('ky');
% colormap('gray');


U_2D_score = zeros(256,256); 
% To find how many times that point along kx, ky appears in the x% of points 
% contributes to the maximum energy of the signal

U_2D_points = zeros(1,172);
% To find how many points are needed actually to contribute to tol% of the
% energy of the signal,can help in determining the sparsity of the solution
% per slice in case of 2D constant undersampling matrix

tol = 0.99; % percent of energy needs to be covered, i.e what kind of image are we
% accepting to be an ideal image. This means that any image which contains
% greater than or equal to 100*tol % of the energy of the original is considered to be as 
% good as the original

for slice_no = 1:172
    slice_no
    temp_img = fImage(:,:,slice_no);
    temp_ft= ft(temp_img);
    temp_ft_vec = temp_ft(:);
    
    
        max_energy_points = findnpoints(temp_ft_vec,tol);
        U_2D_points(1,slice_no) = length(max_energy_points);
        U_2D_score(max_energy_points) = U_2D_score(max_energy_points)+1;
  
end

% U_2D_score includes just tol
U_2D_score = U_2D_score/(max(max(U_2D_score))); % Normalizing the score to find out 
% that particular position occurs how much percent of time,172 for no of
% slices and 256 for rows in each slice

guess_error = 0.05; % We are ok with our asssumptions (look at tol) are wrong

% Find all points that justify our assumptions of error tolerance and
% correctness of guess

U_2D_score_all = zeros(size(U_2D_score));
% including all assumptions of guess_error and tol
pos_u = find(U_2D_score> guess_error);
U_2D_score_all(pos_u) = U_2D_score(pos_u);

min_sparsity = sum(sum(round(U_2D_score_all>0)))/(256*256);


figure(2);
% hold on;
subplot(1,2,1);
imagesc(U_2D_score);
title(['Frequency of occurance of points in kx space, tol = ' num2str(tol)]);
subplot(1,2,2);
imagesc(U_2D_score_all);
colormap('gray');
% hold off;
title(['Frequency of occurance of points in kx space, tol = ' num2str(tol) ' guess\_err = ' num2str(guess_error)]);

% Defining the undersampling pattern U
U = round(U_2D_score_all>0);

figure(3);
imagesc(U);
title(['U for tol = ' num2str(tol) ', guess\_err = ' num2str(guess_error)]);
colormap('gray');

U2DImage = zeros(size(fImage));
% 

err_vec=[];
count_vec = [];

for slice_no=1:172
slice_no
ground_truth = fImage(:,:,slice_no);

% COMPRESSED SENSING

% We will undersample in FT
gt_s = ft(ground_truth);

% undersampled data
U_sampled_data = gt_s.*U;

% CS_reconstruction
imhat = zeros(size(U_sampled_data));
imhat = abs(fftshift(ifft2(U.*gt_s + (1-U).*fft2(imhat))));
imhat_first = imhat; % The undersampled image
err = imrotate(abs(imhat),0)/norm(imhat(:))-ground_truth/norm(ground_truth(:));
err_per = norm(err(:));
count = 0;

while (err_per>1-tol && count<100) % threshold error = 1%, max_iterations = 100

% enforcing sparsity in wavelet domain
wavelet_used = 'db2';
n_levels = 2;
[c s] = wavedec2(imhat,n_levels,wavelet_used);
[val pos] = sort(abs(c));
val(1:round(0.25*length(imhat(:)))) = 0; % Retaining only 75% of values
c_new = zeros(size(c));
for j=1:length(c)
    c_new(pos(j)) = c(pos(j));
end
imhat = waverec2(c,s,wavelet_used);


% data consistency
imhat = abs((ift(U.*gt_s + (1-U).*ft(imhat)))); 

% error 
err = imrotate(abs(imhat),0)/norm(imhat(:))-ground_truth/norm(ground_truth(:));
err_per = norm(err(:));
% iteration update
count = count+1;
end
err_vec = [err_vec;err_per];
count_vec = [count_vec; count];

U2DImage(:,:,slice_no) = abs(imhat);
UnSampimage(:,:,slice_no) = imhat_first;
end

figure(4);
subplot(2,1,1);
plot(err_vec);
title('reconstruction error vs slice no');
subplot(2,1,2);
plot(count_vec);
title('no of iter to converge (max = 100)');

figure(5);
imshow3Dfull(fImage);
title('Original data');
figure(6);
imshow3Dfull(UnSampimage);
title('Undersampled data');
figure(7);
imshow3Dfull(U2DImage);
title('Undersampled CS recon data');
