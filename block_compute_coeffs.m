function coeffs = block_compute_coeffs(domain, tau, refs, mpxli, mpxlj)
%COMPUTE_COEFFS Summary of this function goes here
%   Detailed explanation goes here

[n, m, o, ~] = size(tau);
[~, p_tot, q_tot] = size(refs);
p = length(mpxli);
q = length(mpxlj);

translated = zeros(p*q,o);
coeffs = zeros(n,m,o);
            
for i=1:n
    for j=1:m
        crt = squeeze(domain(i,j,:,:));
        for k=1:o
            i_range = min(max(mpxli - tau(i,j,k,1), 1), p_tot);
            j_range = min(max(mpxlj - tau(i,j,k,2), 1), q_tot);
            crt_trans = refs(k, i_range, j_range);
            
            translated(:,k) = crt_trans(:);
        end
        coeffs(i,j,:) = inv(translated.'*translated + 1e-3*eye(o)) * translated.' * crt(:);
    end
end
end

