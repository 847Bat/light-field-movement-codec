function [sq, s_range, s_sign, cq, c_range, c_sign, muq, mu_range]  = ...
    residual_compression(ground_truth, predicted, nb_components, nb_bits, quantization_mode)
%RESIUDAL_COMPRESSION Summary of this function goes here
% Use 4 bits, and refer to the fig to have the nb of components for each
% bpp. 
% PCA is a plus because it handles well 2D (no need for zigzag), it
% recognizes patterns in the angular space (= recognize patterns in the 
% image space, ...
% Look at the patterns and the coeffs !

% Do it per block ??? -> some blocks have no residuals...

[M, N, O, P] = size(ground_truth);
residuals = ground_truth - predicted;

ref_rows = zeros(O*P, M*N);
pred_rows = zeros(O*P, M*N);
res_rows = zeros(O*P, M*N);
for i=1:O
    for j=1:P
        ref_rows((i-1)*P+j,:) = reshape(ground_truth(:,:,i,j),[M*N 1]);
        res_rows((i-1)*P+j,:) = reshape(residuals(:,:,i,j),[M*N 1]);
        pred_rows((i-1)*P+j,:) = reshape(predicted(:,:,i,j),[M*N 1]);
    end
end

% PCA
[c, s, ~, ~, ~, mu] = pca(res_rows, 'NumComponents', nb_components);

% Quantization
Q_pca = 2^nb_bits;

s_range = max(abs(s(:))) + 1;
c_range = max(abs(c(:))) + 1;
mu_range = max(abs(mu(:))) + 1;
s_sign = s >= 0;
c_sign = c >= 0;

if quantization_mode == 'lin'
    sq=round(s.*Q_pca./s_range);
    cq=round(c.*Q_pca./c_range);
elseif quantization_mode == 'log'
    scale = exp(Q_pca/2 - 1);   % cause bit of sign
    sq= round(log(abs(s.*scale./s_range) + 1));
    cq= round(log(abs(c.*scale./c_range) + 1));
else
    disp("Unrecognize mode of quantization. Possibilities are lin or log.");
    return;
end

muq = round(mu.*2^15./mu_range);

end

