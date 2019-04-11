function coeffs = compute_coeffs_nosub(domain, tau, refs)
%COMPUTE_COEFFS Summary of this function goes here
%   Detailed explanation goes here

[n, m, o, ~] = size(tau);
[p, q, ~] = size(refs);

translated = zeros(p*q,o);
coeffs = zeros(n,m,o);
            
for i=1:n
    for j=1:m
        crt = im2double(squeeze(domain(i,j,:,:)));
        for k=1:o
            
            crt_trans = zeros(p,q);
            crt_trans(max(1,1+tau(i,j,k,1)):min(p,p+tau(i,j,k,1)), max(1,1+tau(i,j,k,2)):min(q,q+tau(i,j,k,2))) = ...
                refs(max(1,1-tau(i,j,k,1)):min(p,p-tau(i,j,k,1)), max(1,1-tau(i,j,k,2)):min(q,q-tau(i,j,k,2)),k);
            
            % Correction for borders
            if tau(i,j,k,1) < 0
                crt_trans(end + floor(tau(i,j,k,1)) + 1:end,:) = repmat(refs(end,:,k),[abs(tau(i,j,k,1)) 1]);
            elseif tau(i,j,k,1) > 0
                crt_trans(1:floor(tau(i,j,k,1)),:) = repmat(refs(1,:,k),[abs(tau(i,j,k,1)) 1]);
            end
            if tau(i,j,k,2) < 0
                crt_trans(:,end + floor(tau(i,j,k,2)) + 1:end) = repmat(refs(:,end,k),[1 abs(tau(i,j,k,2))]);
            elseif tau(i,j,k,2) > 0
                crt_trans(:,1:floor(tau(i,j,k,2))) = repmat(refs(:,1,k),[1 abs(tau(i,j,k,2))]);
            end
            
            translated(:,k) = crt_trans(:);
        end
        coeffs(i,j,:) = inv(translated.'*translated) * translated.' * crt(:);
    end
end
end

