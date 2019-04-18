function [taus, coeffs, predicted] = predictor_coder(domain, refs, blocks)
%PREDICTOR_CODER Summary of this function goes here
%   Detailed explanation goes here

[M, N, ~, ~] = size(domain);
nb_blocks = size(blocks,1);
Q = size(refs, 1);
coeffs = zeros(nb_blocks,M,N,Q);
taus = zeros(nb_blocks,M,N,Q,2, 'int16');
reverseStr = [];

% Compressing refs TODO : HEVC
Q_refs = 2^16;
refs = floor(refs*Q_refs) / Q_refs;

fprintf('\tFinding translation : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);
    
    crt_ref = refs(:,mpxli,mpxlj);
    crt_domain = domain(:,:,mpxli,mpxlj);

    % Find the translations
    tau = block_predictor_coder(crt_domain, crt_ref, 1);
    taus(i_block, :,:,:,:) = tau;
end
fprintf('\t\tDone\n');

predicted = zeros(size(domain));
reverseStr = [];
fprintf('\tComputing coefficients : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);
    
    tau = squeeze(taus(i_block, :,:,:,:));
    crt_ref = refs(:,mpxli,mpxlj);
    crt_domain = domain(:,:,mpxli,mpxlj);
    
    % Compute coeffs
    coeffs(i_block, :,:,:) = block_compute_coeffs(crt_domain, tau, crt_ref);
end
fprintf('\tDone\n');

% Compressing coeffs
Q_coeffs = 2^16;
coeffs_min = min(coeffs(:));
coeffs_range = max(coeffs(:)) - min(coeffs(:));
coeffs = floor((coeffs - coeffs_min)/coeffs_range*Q_coeffs) * coeffs_range / Q_coeffs + coeffs_min;

end

