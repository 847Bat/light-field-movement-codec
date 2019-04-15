function decoded_LF = decoder(filename)
%DECODER Summary of this function goes here
%   Detailed explanation goes here

% Decompress values
disp("Unpacking data");
load(filename);

blocks = reshape(bi2de(bwunpack(blocks_c, prod(blocks_c_p))), blocks_c_p);
refs = reshape(bi2de(bwunpack(refs_c, prod(refs_c_p))), refs_c_p)/2^16;
taus = reshape(bi2de(bwunpack(taus_c, prod(taus_c_p))), taus_c_p) + cast(taus_c_p2, 'double');
coeffs = reshape(bi2de(bwunpack(coeffs_c, prod(coeffs_c_p))), coeffs_c_p)*coeffs_c_p3/(2^16-1) + coeffs_c_p2;
sq = reshape(bi2de(bwunpack(sq_c, prod(sq_c_p2))), sq_c_p2);
s_sign = reshape(bwunpack(sq_c_p, sq_c_p2(1)), sq_c_p2);
s_range = sq_c_p3;
cq = reshape(bi2de(bwunpack(cq_c, prod(cq_c_p2))), cq_c_p2);
c_sign = reshape(bwunpack(cq_c_p, cq_c_p2(1)), cq_c_p2);
c_range = cq_c_p3;
sizes = sizes_c;
muq = bi2de(bwunpack(muq_c, sizes(1)*sizes(2))).' - 2^15;
muq_range = muq_c_p;

% Predictor decoder
disp("Predictor decoder");
predicted = predictor_decoder(refs, taus, coeffs, blocks);

% Residual decompression
disp("Residuals decompression");
res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign,...
 c_range, muq, muq_range, 4, 'log', sizes);

% Reconstruction
decoded_LF = max(min(predicted + res_reconstructed, 1),0);
% decoded_LF = max(min(predicted, 1),0);

disp("Decoding completed.");
end

