function loc = block_matching(block, refs)
%BLOCK_MATCHING Summary of this function goes here
%   Detailed explanation goes here

[m, n, o, p] = size(block);
[k, o2, p2] = size(refs);
i_range = o2-o+1;
j_range = p2-p+1;

refs_2 = zeros(k, i_range*j_range);
refs_ = zeros(k, o*p, i_range*j_range);
block_2 = zeros(m*n, 1);
block_ = zeros(m*n, o*p);

loc = zeros(m, n, k, 2);
for i_ref=1:k
    for i=1:i_range
        for j=1:j_range
            tmp = refs(i_ref, i:i+o-1, j:j+p-1);
            refs_(i_ref, :, (i-1)*j_range + j) = tmp(:);
            refs_2(i_ref, (i-1)*j_range + j) = tmp(:).'*tmp(:);
        end
    end
end

for i=1:m
    for j=1:n
        tmp = block(i,j,:,:);
        block_((i-1)*n + j,:) = tmp(:);
        block_2((i-1)*n + j) = tmp(:).'*tmp(:);
    end
end

for i_ref=1:k  
    error = block_2 + refs_2(i_ref, :) - ...
        2*block_*reshape(refs_(i_ref, :,:), [o*p i_range*j_range]);
    for i=1:m*n
        mi = min(error(i,:));
        t_tot = find(error(i,:)==mi,1, 'first');
        ty = floor((t_tot-1)/j_range) + 1;
        tx = mod(t_tot-1,j_range) + 1;
        loc(floor((i-1)/n) + 1, mod(i-1,n)+1, i_ref, :) = [ty tx];
    end
end
end

