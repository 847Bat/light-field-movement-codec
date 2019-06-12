function predicted = block_predictor_decoder(refs, tau, coeffs, mpxli, mpxlj, my_closest_refs)
%PREDICTOR_DECODER Summary of this function goes here
%   Detailed explanation goes here

[n, ~, ~] = size(tau);
[o, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);
nb_neighbors = size(my_closest_refs,2);

predicted = zeros(n,p,q);

for i=1:n
    for k=1:nb_neighbors
        i_range = min(max(mpxli - tau(i,k,1), 1), p_tot);
        j_range = min(max(mpxlj - tau(i,k,2), 1), q_tot);
        crt_trans = squeeze(refs(my_closest_refs(i,k), i_range, j_range));
        %crt_ref = squeeze(refs(k, mpxli, mpxlj));

        predicted(i,:,:) = squeeze(predicted(i,:,:)) + coeffs(i,k).*crt_trans;
        %predicted(i,:,:) = squeeze(predicted(i,:,:)) + coeffs(i,k+o).*crt_ref;
    end
end
end

