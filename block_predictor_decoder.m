function predicted = block_predictor_decoder(refs, tau, coeffs, mpxli, mpxlj)
%PREDICTOR_DECODER Summary of this function goes here
%   Detailed explanation goes here

[n, ~, ~] = size(tau);
[o, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);
predicted = zeros(n,p,q);

for i=1:n
    for k=1:o
        i_range = min(max(mpxli - tau(i,k,1), 1), p_tot);
        j_range = min(max(mpxlj - tau(i,k,2), 1), q_tot);
        crt_trans = squeeze(refs(k, i_range, j_range));

        predicted(i,:,:) = squeeze(predicted(i,:,:)) + coeffs(i,k).*crt_trans;
    end
end
end

