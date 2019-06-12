function predicted = predictor_decoder(refs, taus, coeffs, blocks, my_closest_refs)
%PREDICTOR_DECODER Summary of this function goes here
%   Detailed explanation goes here
% Reconstruct predicted image

[~, N, ~, ~] = size(taus);
[~, O, P] = size(refs);
nb_blocks = size(blocks, 1);
predicted = zeros(N, O, P);
reverseStr = [];

fprintf('\tDecoding block : ');
for i_block = 1:nb_blocks
    msg = sprintf('%d/%d', i_block, nb_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    mpxli = blocks(i_block,1,1):blocks(i_block,1,2);
    mpxlj = blocks(i_block,2,1):blocks(i_block,2,2);

    crt_tau = reshape(taus(i_block, :,:,:), [size(taus,2), size(taus,3), size(taus,4)]);
    crt_coeffs = reshape(coeffs(i_block, :,:), [size(coeffs,2), size(coeffs,3)]);

    predicted(:,mpxli,mpxlj) = ...
        max(min(block_predictor_decoder(refs, crt_tau, crt_coeffs, mpxli, mpxlj, my_closest_refs),1),0);
end
fprintf('\tDone\n');

end

