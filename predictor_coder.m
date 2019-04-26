function [taus, coeffs] = predictor_coder(domain, refs, blocks)
%PREDICTOR_CODER Summary of this function goes here
%   Detailed explanation goes here

[N, ~, ~] = size(domain);
nb_blocks = size(blocks,1);
[Q, O, P] = size(refs);
coeffs = zeros(nb_blocks,N,Q);
taus = zeros(nb_blocks,N,Q,2);
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
    mpxli_e = max(blocks(i_block,1,1)-10,1):min(blocks(i_block,1,2)+10, O);
    mpxlj_e = max(blocks(i_block,2,1)-10,1):min(blocks(i_block,2,2)+10, P);    
    
    crt_domain = domain(:,mpxli,mpxlj);
    % crt_ref = refs(:,mpxli_e,mpxlj_e);
    crt_ref = refs(:,mpxli,mpxlj);

    % Find the translations
    tau = block_predictor_coder(crt_domain, crt_ref, 1);
%     tau = - (block_matching(crt_domain, crt_ref) - ...
%         reshape([mpxli(1) - mpxli_e(1) mpxlj(1) - mpxlj_e(1)], [1 1 1 2]) - 1);
    taus(i_block, :,:,:) = tau;
end
fprintf('\t\tDone\n');

reverseStr = [];
fprintf('\tComputing coefficients : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    %reverseStr = repmat(sprintf('\b'), 1, length(msg));

    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);
    
    tau = squeeze(taus(i_block, :,:,:));
    crt_domain = domain(:,mpxli,mpxlj);
    
    % Compute coeffs
    crt_coeffs = block_compute_coeffs(crt_domain, tau, refs, mpxli, mpxlj);
    
    coeffs(i_block,:,:) = crt_coeffs;
    
end
fprintf('\tDone\n');


% Compressing coeffs
Q_coeffs = 2^16;
coeffs_min = min(coeffs(:));
coeffs_range = max(coeffs(:)) - min(coeffs(:));
coeffs = floor((coeffs - coeffs_min)/coeffs_range*Q_coeffs) * coeffs_range / Q_coeffs + coeffs_min;

end

