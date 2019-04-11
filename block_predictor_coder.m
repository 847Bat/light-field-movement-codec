function tau = block_predictor_coder(domain, refs, prec)
%PREDICTOR Summary of this function goes here
%   Detailed explanation goes here

[n, m, ~, ~] = size(domain);
[~, ~, o] = size(refs);

f_refs = zeros(size(refs));
for i=1:o
    f_refs(:,:,i) = fft2(refs(:,:,i));
end

tau = zeros(n,m, o, 2);
for i=1:n
    for j=1:m
        crt = squeeze(domain(i,j,:,:));
        f_crt = fft2(crt);
        for k=1:o
            tau(i,j,k,:) = find_translation(f_crt,f_refs(:,:,k),prec);
        end
    end
end

end

