function loc = block_predictor_coder(domain, refs, closest_refs_grey, prec)
%PREDICTOR Summary of this function goes here
%   Detailed explanation goes here

[n, ~, ~] = size(domain);
[o, ~, ~] = size(refs);
nb_neighbors = size(closest_refs_grey,2);

f_refs = zeros(size(refs));
for i=1:o
    f_refs(i,:,:) = fft2(squeeze(refs(i,:,:)));
end

loc = zeros(n, nb_neighbors, 2);
for i=1:n
    crt = squeeze(domain(i,:,:));
    f_crt = fft2(crt);
    for k=1:nb_neighbors
        loc(i,k,:) = find_translation(f_crt,squeeze(f_refs(closest_refs_grey(i,k),:,:)),prec);
    end
end

end

