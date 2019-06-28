function coeffs = block_compute_coeffs(domain, tau, refs, mpxli, mpxlj, closest_refs_grey)
%COMPUTE_COEFFS Summary of this function goes here
%   Detailed explanation goes here

[n, ~, ~] = size(tau);
[~, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);
nb_neighbors = size(closest_refs_grey,2);

translated = zeros(p*q,nb_neighbors);
coeffs = zeros(n,nb_neighbors);
            
for i=1:n
    crt = squeeze(domain(i,:,:));
    for k=1:nb_neighbors
        if floor(tau(i,k,1)) ~= tau(i,k,1) || floor(tau(i,k,2)) ~= tau(i,k,2)
            i_range = min(max(mpxli - floor(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - floor(tau(i,k,2)), 1), q_tot);
            trans_ff = refs(closest_refs_grey(i,k), i_range, j_range);

            i_range = min(max(mpxli - ceil(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - ceil(tau(i,k,2)), 1), q_tot);
            trans_cc = refs(closest_refs_grey(i,k), i_range, j_range);

            i_range = min(max(mpxli - floor(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - ceil(tau(i,k,2)), 1), q_tot);
            trans_fc = refs(closest_refs_grey(i,k), i_range, j_range);

            i_range = min(max(mpxli - ceil(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - floor(tau(i,k,2)), 1), q_tot);
            trans_cf = refs(closest_refs_grey(i,k), i_range, j_range);

            alpha_i = ceil(tau(i,k,1)) - tau(i,k,1);
            alpha_j = ceil(tau(i,k,2)) - tau(i,k,2);
            crt_trans = alpha_j*(alpha_i*trans_ff + (1 - alpha_i)*trans_cf) + ...
                (1 - alpha_j)*(alpha_i*trans_fc + (1 - alpha_i)*trans_cc);
        else
            i_range = min(max(mpxli - floor(tau(i,k,1)), 1), p_tot);
            j_range = min(max(mpxlj - floor(tau(i,k,2)), 1), q_tot);
            crt_trans = refs(closest_refs_grey(i,k), i_range, j_range);
        end

        translated(:,k) = crt_trans(:);
    end
    coeffs(i,:) = inv(translated.'*translated + 1e-3*eye(nb_neighbors)) * translated.' * crt(:);
end
end

