%%
addpath(genpath('modules/subaxis'));
addpath(genpath('modules/light-field-toolbox'));
addpath(genpath('modules/light-field-graph-codec'));

load 'data/4DLF/People/Fountain_&_Vincent_2.mat'

filename = 'compressed.mat';

% Transform to double
LF_double = im2double(LF);
grey_LF = 0.299 * LF_double(:,:,:,:,1) + 0.587 * LF_double(:,:,:,:,2) + 0.114 * LF_double(:,:,:,:,3);
grey_LF = floor(grey_LF*(2^16))./2^16;


%% Encoding
encoder(grey_LF, filename);

%% Display efficiency
file = dir(filename);
disp("Final bpp : " + num2str(file.bytes*8/numel(grey_LF)));
stats = whos('-file',filename);
bpp = num2cell([stats.bytes].'*8/numel(grey_LF));
[stats.bpp] = bpp{:};

for i = 1:length(stats)
    if stats(i).bpp > 1e-5
        msg1 = sprintf('%s',stats(i).name);
        tab = repmat(sprintf('\t'), 1, 2 - fix(length(msg1)/8));
        msg2 =sprintf('%.5f\n',stats(i).bpp);
        fprintf([msg1, tab, msg2]);
    end
end

%% Decoding
decoded_LF = decoder(filename);

%% Assess quality
disp("mean psnr :");
disp(lf_psnr(decoded_LF, grey_LF));

%%
LFDispVidCirc(repmat(decoded_LF, [1 1 1 1 3]));
