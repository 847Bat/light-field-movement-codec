function loc = block_predictor_coder(domain, refs, prec)
%PREDICTOR Summary of this function goes here
%   Detailed explanation goes here

[n, m, ~, ~] = size(domain);
[o, p_tot, q_tot] = size(refs);

f_refs = zeros(size(refs));
for i=1:o
    f_refs(i,:,:) = fft2(squeeze(refs(i,:,:)));
end

loc = zeros(n,m, o, 2);
for i=1:n
    for j=1:m
        crt = squeeze(domain(i,j,:,:));
        f_crt = fft2(crt, size(refs,2), size(refs,3));
        for k=1:o
%             R = squeeze(f_refs(k,:,:));
%             T = f_crt;
%             cross = abs(ifft2(R.*conj(T)./abs(R.*conj(T))));
%             
%             ma = max(cross(:));
%             [t1, t2] = find(cross==ma,1,'first');
%             loc(i,j,k,:) = [t1 t2];
            
            loc(i,j,k,:) = find_translation(f_crt,squeeze(f_refs(k,:,:)),prec);
        end
    end
end

end

