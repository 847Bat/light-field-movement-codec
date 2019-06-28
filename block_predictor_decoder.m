function predicted = block_predictor_decoder(refs, tau, coeffs, mpxli, mpxlj, my_closest_refs)
%PREDICTOR_DECODER Summary of this function goes here
%   Detailed explanation goes here

[n, ~, ~] = size(tau);
[~, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);
nb_neighbors = size(my_closest_refs,2);

predicted = zeros(n,p,q);

for i=1:n
    for k=1:nb_neighbors
        if floor(tau(i,k,1)) ~= tau(i,k,1) || floor(tau(i,k,2)) ~= tau(i,k,2)
            i_range = min(max(mpxli - floor(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - floor(tau(i,k,2)), 1), q_tot);
            trans_ff = squeeze(refs(my_closest_refs(i,k), i_range, j_range));

            i_range = min(max(mpxli - ceil(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - ceil(tau(i,k,2)), 1), q_tot);
            trans_cc = squeeze(refs(my_closest_refs(i,k), i_range, j_range));

            i_range = min(max(mpxli - floor(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - ceil(tau(i,k,2)), 1), q_tot);
            trans_fc = squeeze(refs(my_closest_refs(i,k), i_range, j_range));

            i_range = min(max(mpxli - ceil(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - floor(tau(i,k,2)), 1), q_tot);
            trans_cf = squeeze(refs(my_closest_refs(i,k), i_range, j_range));

            alpha_i = ceil(tau(i,k,1)) - tau(i,k,1);
            alpha_j = ceil(tau(i,k,2)) - tau(i,k,2);
            crt_trans = alpha_j*(alpha_i*trans_ff + (1 - alpha_i)*trans_cf) + ...
                (1 - alpha_j)*(alpha_i*trans_fc + (1 - alpha_i)*trans_cc);
        else
            i_range = min(max(mpxli - floor(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - floor(tau(i,k,2)), 1), q_tot);
            crt_trans = squeeze(refs(my_closest_refs(i,k), i_range, j_range));
        end

        predicted(i,:,:) = squeeze(predicted(i,:,:)) + coeffs(i,k).*crt_trans;
    end
end
end

