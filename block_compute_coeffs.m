function coeffs = block_compute_coeffs(domain, tau, refs)
%COMPUTE_COEFFS Summary of this function goes here
%   Detailed explanation goes here

[n, m, o, ~] = size(tau);
[~, p, q] = size(refs);

translated = zeros(p*q,o);
coeffs = zeros(n,m,o);
            
for i=1:n
    for j=1:m
        crt = squeeze(domain(i,j,:,:));
        for k=1:o
            
            crt_trans = zeros(p,q);
            crt_trans(max(1,1+tau(i,j,k,1)):min(p,p+tau(i,j,k,1)), max(1,1+tau(i,j,k,2)):min(q,q+tau(i,j,k,2))) = ...
                refs(k,max(1,1-tau(i,j,k,1)):min(p,p-tau(i,j,k,1)), max(1,1-tau(i,j,k,2)):min(q,q-tau(i,j,k,2)));
            
            % Correction for borders
            if tau(i,j,k,1) < 0
                crt_trans(end + floor(tau(i,j,k,1)) + 1:end,:) = repmat(refs(k,end,:),[1 abs(tau(i,j,k,1)) 1]);
            elseif tau(i,j,k,1) > 0
                crt_trans(1:floor(tau(i,j,k,1)),:) = repmat(refs(k,1,:),[1 abs(tau(i,j,k,1)) 1]);
            end
            if tau(i,j,k,2) < 0
                crt_trans(:,end + floor(tau(i,j,k,2)) + 1:end) = repmat(refs(k,:,end),[1 1 abs(tau(i,j,k,2))]);
            elseif tau(i,j,k,2) > 0
                crt_trans(:,1:floor(tau(i,j,k,2))) = repmat(refs(k,:,1),[1 1 abs(tau(i,j,k,2))]);
            end
            
            translated(:,k) = crt_trans(:);
        end
        coeffs(i,j,:) = inv(translated.'*translated + 1e-3*eye(o)) * translated.' * crt(:);
    end
end
end

