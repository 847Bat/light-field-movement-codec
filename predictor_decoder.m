function predicted = predictor_decoder(refs, taus, coeffs, blocks)
%PREDICTOR_DECODER Summary of this function goes here
%   Detailed explanation goes here
% Reconstruct predicted image

[~, M, N, ~, ~] = size(taus);
[~, O, P] = size(refs);
nb_blocks = size(blocks, 1);
predicted = zeros(M, N, O, P);
reverseStr = [];

fprintf('\tDecoding block : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);
    
    crt_ref = refs(:,mpxli,mpxlj);

    crt_tau = squeeze(taus(i_block, :,:,:,:));
    crt_coeffs = squeeze(coeffs(i_block, :,:,:));

    predicted(:,:,mpxli,mpxlj) = ...
        max(min(block_predictor_decoder(crt_ref, crt_tau, crt_coeffs),1),0);
end
fprintf('\tDone\n');

end

