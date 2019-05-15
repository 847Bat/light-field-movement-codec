%%
% main10 profile
% ffmpeg ppm to yuv444 420 or 422
% hm software hhi hm-16.9

%%
addpath(genpath('modules/subaxis'));
addpath(genpath('modules/light-field-toolbox'));
addpath(genpath('modules/light-field-graph-codec'));
addpath(genpath('modules/03_Scripts'));

img = 'I09';

filename = [img 'compressed.mat'];

%% Build masks and options

% Attention : use non-translated refs ? -> no

LF = load_LF(img);
LF = LF(2:14, 2:14, :, :, :);
LF_10b = uint16(double(LF)/65535*1023);

% 420 10 [15] + 1:4, +, and 2:10  -> 0.74/Y=42.172 !!! (nt yes)

% 422 12 [15] + idem -> 0.71/YUV=40.03 !!! (nt no)   
% 422 15 [20] + idem
% 444 15 [10] + idem -> 0.72/YUV=39.85

mask_refs = zeros(size(LF, 1), size(LF, 2), 'logical');
mask_refs(1:4:end,1:4:end) = 1;
mask_refs(7,1:6:end) = 1;
mask_refs(1:6:end,7) = 1;
% mask_refs(1:12:end,1:12:end) = 0;
mask_refs(2:10:end,2:10:end) = 1;

masks = zeros(1, size(LF,1), size(LF,2), 'logical');
masks(1, :, :) = 1;
masks(1, mask_refs) = 0;

nb_components = [10];

qp = 15;        % Transparent until 18/24 ?
fmt = '444';

%% Encoding
status = encoder(LF_10b, filename, mask_refs, masks, nb_components, qp, fmt, img);

%% Display efficiency

mask1 = squeeze(mask_refs(:,:,1,1));
mask_4 = zeros(4, fix(numel(mask1)/4 + 1));
mask_4(1:numel(mask1)) = mask1(:);

mask_name = '';
for i=1:size(mask_4,2)
    mask_name = [mask_name dec2hex(bin2dec(int2str(mask_4(:,i)).'))];
end
    
global_folder = ['refs/' img '/' mask_name];
folder = ['refs/' img '/' mask_name '/' fmt '/' int2str(qp)];

file = dir(filename);
refs_file = dir([folder '/refs.bin']);
disp("Final bpp : " + num2str((file.bytes + refs_file.bytes)*8/(13*13*434*625)));
disp("bpp of refs : " + num2str(refs_file.bytes*8/(13*13*434*625)));
stats = whos('-file',filename);
bpp = num2cell([stats.bytes].'*8/(13*13*434*625));
[stats.bpp] = bpp{:};

for i = 1:length(stats)
    if stats(i).bpp > 1e-4
        msg1 = sprintf('%s',stats(i).name);
        tab = repmat(sprintf('\t'), 1, 2 - fix(length(msg1)/8));
        msg2 =sprintf('%.5f\n',stats(i).bpp);
        fprintf([msg1, tab, msg2]);
    end
end

%% Decoding
decoded_LF = decoder(LF_10b, filename);

%% Assess quality
disp("PSNR :");

[PSNR_Y, PSNR_U, PSNR_V, PSNR_YUV, PSNR_Y_mean, PSNR_U_mean, PSNR_V_mean, PSNR_YUV_mean] = ComputePSNR(decoded_LF, LF_10b)

%%
LFDispVidCirc(im2double(decoded_LF));
