function res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign, c_range, muq, mu_range, n_bits, quantization_mode, sizes)
%RESIDUAL_DECOMPRESSION Summary of this function goes here
%   Detailed explanation goes here
M = sizes(1);
N = sizes(2);
O = sizes(3);
P = sizes(4);
Q_pca = 2^n_bits;

sq = sq + .5;
cq = cq + .5;

mu = muq * mu_range/2^15;
if quantization_mode == 'lin'
    s = (s_sign*2-1) .* sq*2/Q_pca*s_range;
    c = (c_sign*2-1) .* cq*2/Q_pca*c_range;
elseif quantization_mode == 'log'
    scale = exp(Q_pca/2 - 1);   % cause bit of sign
    s = (s_sign*2-1) .* (exp(sq) - 1)/scale * s_range;
    c = (c_sign*2-1) .* (exp(cq) - 1)/scale * c_range;
else
    disp("Unrecognize mode of quantization. Possibilities are lin or log.");
    return;
end

res_rows_reconstructed = s*c.' + repmat(mu,length(s),1);
res_reconstructed = zeros(M,N,O,P);
for i=1:O
    for j=1:P
        res_reconstructed(:,:,i,j) = reshape(res_rows_reconstructed((i-1)*P+j,:), [M N]);
    end
end
end

