%% PCA reconstruction / compressions
clear all;
close all;
clc;


% INPUT
for i=1:172
    str=strcat('s6192/','*.MRDC.',num2str(i));
    listing(i) = dir(str);
end

fImage=zeros(256,256,length(listing));
for i=1:length(listing)
    fImage(:,:,i)=dicomread(strcat('s6192/',listing(i).name));
end

block_size = 16;
% partitioning the ft domain into blocks of size 16x16. Our image is
% 256x256 so the number of blocks would be 256x256/16x16 = 256;
n_blocks = 256*256/(block_size)^2;


thresh = 0.9; % retaining the amount of data

U_cell = cell(1,n_blocks);
S_cell = cell(1,n_blocks);
V_cell = cell(1,n_blocks);
orig_data = cell(1,n_blocks);
image_test = zeros(256,256);
k_vec = zeros(1, n_blocks); % no of basis vectors to consider for each block
energy_vec = zeros(1, n_blocks);

for slice_no = 1:172
    slice_no
    temp_img = fImage(:,:,slice_no);
    temp_ft= ((temp_img));
    
    for block_no = 1:n_blocks
        top_left_x = 16*(ceil(block_no/16)-1)+1;
        top_left_y = mod(16*(block_no-1)+1,256);
        curr_img = temp_ft(top_left_x:top_left_x+15, top_left_y:top_left_y+15);
        curr_img = curr_img(:);
        temp_data = orig_data{block_no};
        temp_data = [temp_data curr_img];
        orig_data{block_no} = temp_data;
    end
    
end

for block_no = 1:n_blocks
    curr_block_data = orig_data{block_no};
    [U, S, V] = svd(curr_block_data);
    U_cell{block_no} = U;
    S_cell{block_no} = S;
    V_cell{block_no} = V;
    energy_vec(block_no) = norm(orig_data{block_no},'fro');
end

energy_vec = energy_vec/ sum(energy_vec);

override_variable = zeros(1, n_blocks); % To override and sample the whole block

[val pos]  = sort(energy_vec, 'descend');

override_variable(pos(1)) = 256; % To sample the whole blocks with the highest energy (we are sampling the highest 16 here)

imag_test = zeros(16,16);
energy_block = zeros(256, 256);
U_sampling_pattern = zeros(256, 256);

for block_no = 1:n_blocks
        
    diag_s = diag(S_cell{block_no});
    diag_cdf = my_cdf(diag_s);
    [pos] = find(diag_cdf >= thresh);
    if length(pos) > 0
    k_vec(block_no) = pos(1);
    end
    energy_block(:, block_no) = sum(abs(orig_data{block_no}),2);
    
    curr_energy_vec = energy_block(:, block_no);
    
    block_sample_pattern = zeros(size(curr_energy_vec));
    [val pos] = sort(curr_energy_vec,'descend');
    block_sample_pattern(pos(1:(max(k_vec(block_no),override_variable(block_no))))) = 1;
    block_sample_pattern = reshape(block_sample_pattern, [16 16]);
    
    top_left_x = 16*(ceil(block_no/16)-1)+1;
    top_left_y = mod(16*(block_no-1)+1,256);
    U_sampling_pattern(top_left_x:top_left_x+15, top_left_y:top_left_y+15) = block_sample_pattern;
    
    
    
    
%     figure(ceil(block_no/64));
%     subplot(8,8,mod(block_no-1,64)+1);
%     imagesc(reshape(energy_block(:, block_no), [16 16]));
%     
    
end


figure; plot(k_vec);
title('no of basis vectors (k\_j) vs block B\_j');

figure; imagesc(reshape(energy_vec, [16, 16]));
title('Avg energy distribution across Blocks');



%% Reconstruction using PCA basis

fImage_pca_recon = zeros(size(fImage));
err_vec = [];

for slice_no = 1:172
    
ground_truth = fImage(:,:,slice_no);
gt_s = (ground_truth);
imhat_first = zeros(256, 256);

for block_no = 1:n_blocks
        top_left_x = 16*(ceil(block_no/16)-1)+1;
        top_left_y = mod(16*(block_no-1)+1,256);
        curr_img = gt_s(top_left_x:top_left_x+15, top_left_y:top_left_y+15);
        curr_est = image_est(curr_img, U_cell{block_no}, k_vec(block_no));
        imhat_first(top_left_x:top_left_x+15, top_left_y:top_left_y+15) = curr_est;
end

imhat_first = abs((imhat_first));

fImage_pca_recon(:,:,slice_no) = imhat_first;
err = imrotate(abs(imhat_first),0)/norm(imhat_first(:))-ground_truth/norm(ground_truth(:));
err_per = norm(err(:));
err_vec = [err_vec err_per];

end

figure(5);
imshow3Dfull(fImage);
title('Original data');
figure(6);
imshow3Dfull(fImage_pca_recon);
title('Undersampled data');


figure;
plot(err_vec);
title(['error vs slice no at, ' '% data retained = ' num2str(100*thresh) ', sparsity = ' num2str(sum(k_vec)/172/256)]);