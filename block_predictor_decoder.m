function predicted = block_predictor_decoder(refs, tau, coeffs, mpxli, mpxlj)
%PREDICTOR_DECODER Summary of this function goes here
%   Detailed explanation goes here

[n, m, ~, ~] = size(tau);
[o, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);
predicted = zeros(n,m,p,q);

for i=1:n
    for j=1:m
        for k=1:o
            i_range = min(max(mpxli - tau(i,j,k,1), 1), p_tot);
            j_range = min(max(mpxlj - tau(i,j,k,2), 1), q_tot);
            crt_trans = squeeze(refs(k, i_range, j_range));
            
            predicted(i,j,:,:) = squeeze(predicted(i,j,:,:)) + coeffs(i,j,k).*crt_trans;
        end
    end
end
end

