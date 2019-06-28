%% Add paths of functions
addpath(genpath('modules/light-field-toolbox'));
addpath(genpath('modules/03_Scripts'));

img = 'I01';

filename = [img 'compressed'];


LF = load_LF(img, 15, 15, 434, 626);
LF = LF(2:14, 2:14, :, :, :);
LF_10b = uint16(double(LF)/65535*1023);

%% Build masks and options (configuration 2, 0.75bpp)

mask_refs = zeros(size(LF, 1), size(LF, 2), 'logical');
mask_refs(1:2:end,1:2:end) = 1;
mask_refs(2:2:end,2:2:end) = 1;

masks = zeros(1, size(LF,1), size(LF,2), 'logical');
masks(1, :, :) = 1;
masks(1, mask_refs) = 0;

nb_neighbors = [[20 4 4]];
nb_components = [0];
qp = 25;
fmt = '420';
nb_blocks = 70;

%% Build masks and options (configuration 1, 0.75bpp)

% mask_refs = zeros(size(LF, 1), size(LF, 2), 'logical');
% mask_refs(1:4:end,1:4:end) = 1;
% mask_refs(7,1:6:end) = 1;
% mask_refs(1:6:end,7) = 1;
% mask_refs(2:10:end,2:10:end) = 1;
% 
% masks = zeros(1, size(LF,1), size(LF,2), 'logical');
% masks(1, :, :) = 1;
% masks(1, mask_refs) = 0;
% 
% nb_neighbors = [[25 4 4]];
% nb_components = [20];
% qp = 12;
% fmt = '420';
% nb_blocks = 100;

%% Encoding
status = encoder(LF_10b, filename, mask_refs, masks, nb_components, qp, ...
        fmt, img, nb_neighbors, nb_blocks);

%% Display bitrates

mask1 = squeeze(mask_refs(:,:,1,1));
mask_4 = zeros(4, fix(numel(mask1)/4 + 1));
mask_4(1:numel(mask1)) = mask1(:);

mask_name = '';
for j=1:size(mask_4,2)
    mask_name = [mask_name dec2hex(bin2dec(int2str(mask_4(:,j)).'))];
end
    
global_folder = ['refs/' img '/' mask_name];
folder = ['refs/' img '/' mask_name '/' fmt '/' int2str(qp)];

file_coeffs = dir([filename '_coeffs.mat']);
file_params = dir([filename '_params.mat']);
file_residuals = dir([filename '_residuals.mat']);
refs_file = dir([folder '/refs.bin']);
disp("Final bpp : " + num2str((file_coeffs.bytes + file_params.bytes + file_residuals.bytes+ refs_file.bytes)*8/(13*13*434*625)));
disp("bpp of refs : " + num2str(refs_file.bytes*8/(13*13*434*625)));
disp("bpp of coeffs : " + num2str(file_coeffs.bytes*8/(13*13*434*625)));
disp("bpp of residuals : " + num2str(file_residuals.bytes*8/(13*13*434*625)));
disp("bpp of params : " + num2str(file_params.bytes*8/(13*13*434*625)));


%% Decoding
decoded_LF = decoder(filename);

%% Assess quality
disp("Quality metrics :");

[PSNR_Y, PSNR_U, PSNR_V, PSNR_YUV, PSNR_Y_mean, PSNR_U_mean, PSNR_V_mean, PSNR_YUV_mean] = ComputePSNR(decoded_LF, LF_10b);
%[SSIM_Y, SSIM_U, SSIM_V, SSIM_YUV, SSIM_Y_mean, SSIM_U_mean, SSIM_V_mean, SSIM_YUV_mean] = ComputeSSIM(decoded_LF, LF_10b);
PSNR_Y_mean
PSNR_YUV_mean
%SSIM_Y_mean

%%
LFDispVidCirc(im2double(decoded_LF));
