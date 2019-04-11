function [taus, coeffs] = predictor_coder(grey_LF, refs, blocks)
%PREDICTOR_CODER Summary of this function goes here
%   Detailed explanation goes here

[M, N, ~, ~] = size(grey_LF);
nb_blocks = size(blocks,1);
Q = size(refs, 3);
coeffs = zeros(nb_blocks,M,N,Q);
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
    tau = block_predictor_coder(domain, crt_ref, 1);

    % Compute coeffs
    coeffs(i_block, :,:,:) = block_compute_coeffs(domain, tau, crt_ref);
    taus(i_block, :,:,:,:) = tau;
end
fprintf('\tDone\n');

end

