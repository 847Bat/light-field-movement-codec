function decoded_LF = decoder(filename)
%DECODER Summary of this function goes here
%   Detailed explanation goes here

% Decompress values
disp("Unpacking data");
load(filename);

blocks = reshape(bi2de(bwunpack(blocks_c, prod(blocks_c_p))), blocks_c_p);
refs = reshape(bi2de(bwunpack(refs_c, prod(refs_c_p))), refs_c_p)/2^16;
masks = reshape(bwunpack(masks_c, prod(masks_c_p)), masks_c_p);
mask_refs = reshape(bwunpack(mask_refs_c, prod(mask_refs_c_p)), mask_refs_c_p);

decoded_LF = zeros(size(mask_refs));
decoded_LF(mask_refs) = refs(:);

for i=1:size(masks,1)
    disp("LEVEL " + int2str(i));
    
    taus = reshape(bi2de(bwunpack(taus_c{i}, prod(taus_c_p{i}))), taus_c_p{i}) + cast(taus_c_p2{i}, 'double');
    coeffs = reshape(bi2de(bwunpack(coeffs_c{i}, prod(coeffs_c_p{i}))), coeffs_c_p{i})*coeffs_c_p3{i}/(2^16-1) + coeffs_c_p2{i};
    sq = reshape(bi2de(bwunpack(sq_c{i}, prod(sq_c_p2{i}))), sq_c_p2{i});
    s_sign = reshape(bwunpack(sq_c_p{i}, sq_c_p2{i}(1)), sq_c_p2{i});
    s_range = sq_c_p3{i};
    cq = reshape(bi2de(bwunpack(cq_c{i}, prod(cq_c_p2{i}))), cq_c_p2{i});
    c_sign = reshape(bwunpack(cq_c_p{i}, cq_c_p2{i}(1)), cq_c_p2{i});
    c_range = cq_c_p3{i};
    sizes = sizes_c{i};
    muq = bi2de(bwunpack(muq_c{i}, sizes(1))).' - 2^15;
    muq_range = muq_c_p{i};
    
    % Predictor decoder
    disp("Predictor decoder");
    predicted = predictor_decoder(refs, taus, coeffs, blocks);
    
    % Residual decompression
    disp("Residuals decompression");
    res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign,...
     c_range, muq, muq_range, 4, 'log', sizes);
 
    % Getting ready for another turn
    refs_obtained = max(min(predicted + res_reconstructed, 1),0);
    refs = cat(1, refs_obtained, refs);
    
    % Stock that
    decoded_LF(masks(i,:,:,:,:)) = refs_obtained(:);
end


disp("Decoding completed.");
end

