function coeffs = block_compute_coeffs(domain, tau, refs, mpxli, mpxlj, closest_refs_grey)
%COMPUTE_COEFFS Summary of this function goes here
%   Detailed explanation goes here

[n, o, ~] = size(tau);
[~, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);
nb_neighbors = size(closest_refs_grey,2);

translated = zeros(p*q,nb_neighbors);%2*o);
coeffs = zeros(n,nb_neighbors);%2*o);
            
for i=1:n
    crt = squeeze(domain(i,:,:));
    for k=1:nb_neighbors
        i_range = min(max(mpxli - tau(i,k,1), 1), p_tot);
        j_range = min(max(mpxlj - tau(i,k,2), 1), q_tot);
        crt_trans = refs(closest_refs_grey(i,k), i_range, j_range);
        %crt_ref = refs(closest_refs_grey(i,k, mpxli, mpxlj);

        translated(:,k) = crt_trans(:);
        %translated(:,k+o) = crt_ref(:);
    end
    coeffs(i,:) = inv(translated.'*translated + 1e-3*eye(o)) * translated.' * crt(:);
end
end

