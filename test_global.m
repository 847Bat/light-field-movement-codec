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

%% Segmentation
resized = imresize(abs(squeeze(grey_LF(8,3,:,:)) - squeeze(grey_LF(8,13,:,:))) + ...
    abs(squeeze(grey_LF(3,8,:,:)) - squeeze(grey_LF(13,8,:,:))), [512 512]);
% resized = imresize(squeeze(grey_LF(8,8,:,:)), [512 512]);
S = qtdecomp(resized,.5, [32 32]);
blocks = repmat(uint8(0),size(S));

for dim = [512 256 128 64 32 16 8 4 2 1]    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

imshow(blocks,[])

%%
my_blocks = zeros(length(find(S>0)), 2, 2);
resizer = [size(grey_LF,3)/512 0 ; 0 size(grey_LF,4)/512];
counter = 1;
for i=1:512
    for j=1:512
        if S(i,j) > 0
            block_512 = [i i + S(i,j) - 1; j j + S(i,j) - 1];
            if ceil((j + S(i,j))*resizer(2,2)) - ceil((j + S(i,j) - 1)*resizer(2,2)) > 1 && j + S(i,j) < 512
                block_512(2,2) = j + S(i,j);
            end
            block_resized = ceil(resizer*block_512);
            my_blocks(counter, :, :) = block_resized;
            counter = counter + 1;
        end
    end
end

%% Encoding
encoder(grey_LF, filename, my_blocks);

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


%% Display coeffs
figure;
m = M;
n = N;
for i=1:m
    for j=1:n
        subaxis(m, n, j, i, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0);
        imagesc(reshape(coeffs(:,i,j,1), [O P])); 
        axis off;
    end
end

figure;
for i_block = 1:nb_blocks
    subaxis(O, P, i_block, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0);
    imagesc(squeeze(coeffs(i_block,:,:,1))); 
    axis off;
end

% -> dct along angular space

%% Display taus

figure;
m = M;
n = N;
for i=1:m
    for j=1:n
        subaxis(m, n, j, i, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0);
        imagesc(reshape(taus(:,i,j,2,2), [P O]).');
        axis off;
    end
end

figure;
m = O;
n = P;
for i_block = 1:nb_blocks
    subaxis(O, P, i_block, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0);
    imagesc(squeeze(taus(i_block,:,:,2,2))); 
    axis off;
end

% -> do clustering along picture space (Warning : need to normalize fft if
% size of blocks is not constant !) and dct along angular space (before
% computing coeffs) (or use coeffs to discard corners ?)

%% Display residuals
residuals = grey_LF - predicted;
range = [min(residuals(:)) max(residuals(:))];

figure;
b=143;
for i=blocks(b,1,1):blocks(b,1,2)
    for j=blocks(b,2,1):blocks(b,2,2)
        subaxis(blocks(b,1,2)-blocks(b,1,1), blocks(b,2,2)-blocks(b,2,1), j-blocks(b,2,1)+1, i-blocks(b,1,1)+1, 'sh', 0.001, 'sv', 0.001, 'padding', 0, 'margin', 0);
        imagesc(squeeze(residuals(:,:,i,j))); 
        axis off;
    end
end

figure;
m = M;
n = N;
for i=1:m
    for j=1:n
        subaxis(m, n, j, i, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0);
        imagesc(squeeze(log(abs(residuals(i,j,:,:))+1))); 
        axis off;
    end
end

figure;
histogram(residuals(:)*65535);  % -> entropy coding, quantization, not so sparse... 
% PCA ? (patterns in space domain so sample is 1 view, features are pxl...)