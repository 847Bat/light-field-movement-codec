%%
addpath(genpath('modules/subaxis'));
addpath(genpath('modules/light-field-toolbox'));
addpath(genpath('modules/light-field-graph-codec'));

load 'data/4DLF/People/Fountain_&_Vincent_2.mat'

%% Object segmentation  (blocks)

LF_double = im2double(LF);
grey_LF = 0.299 * LF_double(:,:,:,:,1) + 0.587 * LF_double(:,:,:,:,2) + 0.114 * LF_double(:,:,:,:,3);

MPXL_size = [31 25];
[M, N, O, P] = size(squeeze(grey_LF));
O = O/MPXL_size(1);
P = P/MPXL_size(2);
nb_blocks = 0*P;
blocks = zeros(nb_blocks, 2, 2);

for i=1:O
    for j=1:P
        blocks((i-1)*P+j,:,:) = [[(i-1)*MPXL_size(1)+1 i*MPXL_size(1)]; [(j-1)*MPXL_size(2)+1 j*MPXL_size(2)]];
    end
end

%% Predictor

% profile on

% Define reference
coded_t = squeeze(grey_LF(3,8,:,:));
coded_r = squeeze(grey_LF(8,13,:,:));
coded_b = squeeze(grey_LF(13,8,:,:));
coded_l = squeeze(grey_LF(8,3,:,:));
coded_m = squeeze(grey_LF(8,8,:,:));

refs = cat(3,coded_t, coded_r, coded_b, coded_l, coded_m);

[taus, coeffs] = predictor_coder(grey_LF, refs, blocks);

% profile off
% profsave

%% Quantization

taus_saved = taus;
coeffs_saved = coeffs;
refs_saved = refs;

Q_taus = 32;
Q_coeffs = 2^8;
Q_refs = 2^16;

taus = min(max(taus_saved,-Q_taus/2),Q_taus/2);   % taus are int, force range

% the rest are float between 0 and 1
coeffs = floor(coeffs_saved*Q_coeffs) / Q_coeffs;
refs = floor(refs_saved*Q_refs) / Q_refs;

%% Decoder

predicted = predictor_decoder(refs, taus, coeffs, blocks);

%% Improvment from linear
mov_err = zeros(M,N);

for i=1:M
    for j=1:N
        ref = squeeze(grey_LF(i+(15-M)/2,j+(15-N)/2,:,:));
        mov_err(i,j) = psnr(squeeze(predicted(i,j,:,:)), ref);
    end
end

%% Efficiency

compress_ratio = 16*size(grey_LF(:)) / (size(taus(:))*log2(Q_taus) + size(coeffs(:))*log2(Q_coeffs) + size(refs(:))*log2(Q_refs))

err_considered = mov_err(2:end-1,2:end-1);
mov_perf = mean(err_considered(:))

%%
LFDispVidCirc(repmat(predicted, [1 1 1 1 3]));

%%
LFDispVidCirc(repmat(lin_pred, [1 1 1 1 3]));

%%
LFDispVidCirc(repmat(grey_LF, [1 1 1 1 3]));

%%

LFDispMousePan(repmat(predicted, [1 1 1 1 3]));

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

