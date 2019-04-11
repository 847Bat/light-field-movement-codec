%%
addpath(genpath('subaxis'));
addpath(genpath('light-field-toolbox'));
addpath(genpath('light-field-graph-codec'));

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

Q = size(refs, 3);
coeffs = zeros(nb_blocks,M,N,Q);
coeffs_lin = zeros(nb_blocks,M,N,Q);
taus = zeros(nb_blocks,M,N,Q,2);
reverseStr = [];

fprintf('\nCompressing block : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);
    
    crt_ref = refs(mpxli,mpxlj,:);
    domain = grey_LF(:,:,mpxli,mpxlj);

    % Find the translations
    tau = predictor_coder(domain, crt_ref, 1);

    % Compute coeffs
    coeffs(i_block, :,:,:) = compute_coeffs_nosub(domain, tau, crt_ref);
    coeffs_lin(i_block, :,:,:) = compute_coeffs_nosub(domain, zeros(size(tau)), crt_ref);
    taus(i_block, :,:,:,:) = tau;
end
fprintf('\tDone\n');

% profile off
% profsave

%% Saving stuff

taus_saved = taus;
coeffs_saved = coeffs;
coeffs_lin_saved = coeffs_lin;
refs_saved = refs;

%% Quantization
Q_taus = 32;
Q_coeffs = 2^8;
Q_refs = 2^16;

taus = min(max(taus_saved,-Q_taus/2),Q_taus/2);   % taus are int, force range

% the rest are float between 0 and 1
coeffs = floor(coeffs_saved*Q_coeffs) / Q_coeffs;
coeffs_lin = floor(coeffs_lin_saved*Q_coeffs) / Q_coeffs;
refs = floor(refs_saved*Q_refs) / Q_refs;

%% Decoder
% Reconstruct predicted image
predicted = zeros(size(grey_LF));
lin_pred = zeros(size(grey_LF));
reverseStr = [];

fprintf('\nDecoding block : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);
    
    crt_ref = refs(mpxli,mpxlj,:);

    crt_tau = squeeze(taus(i_block, :,:,:,:));
    crt_coeffs = squeeze(coeffs(i_block, :,:,:));
    crt_coeffs_lin = squeeze(coeffs_lin(i_block, :,:,:));

    predicted(:,:,mpxli,mpxlj) = ...
        min(predictor_decoder_nosub(crt_ref, crt_tau, crt_coeffs),1);
    lin_pred(:,:,mpxli,mpxlj) = ...
        min(predictor_decoder_nosub(crt_ref, zeros(size(crt_tau)), crt_coeffs_lin),1);
end
fprintf('\tDone\n');

%% Improvment from linear
mov_err = zeros(M,N);
lin_err = zeros(M,N);

for i=1:M
    for j=1:N
        ref = squeeze(grey_LF(i+(15-M)/2,j+(15-N)/2,:,:));
        mov_err(i,j) = psnr(squeeze(predicted(i,j,:,:)), ref);
        lin_err(i,j) = psnr(squeeze(lin_pred(i,j,:,:)), ref);
    end
end

mov_err - lin_err

%% Efficiency

compress_ratio = 16*size(grey_LF(:)) / (size(taus(:))*log2(Q_taus) + size(coeffs(:))*log2(Q_coeffs) + size(refs(:))*log2(Q_refs))

err_considered = mov_err(2:end-1,2:end-1);
mov_perf = mean(err_considered(:))

err_considered = lin_err(2:end-1,2:end-1);
lin_perf = mean(err_considered(:))

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
        subaxis(blocks(b,1,2)-blocks(b,1,1), blocks(b,2,2)-blocks(b,2,1), j-blocks(b,2,1)+1, i-blocks(b,1,1)+1, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0);
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

