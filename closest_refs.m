function my_closest_refs = closest_refs(mask_refs, masks, k)
%CLOSEST_REFS Summary of this function goes here
%   Detailed explanation goes here
crt_refs = mask_refs;
my_closest_refs = {};

for i=1:size(masks,1)
    pos_crt_refs = zeros(sum(sum(crt_refs)), 2);
    pos_crt_domain = zeros(sum(sum(squeeze(masks(i,:,:)))), 2);
    crt_domain = squeeze(masks(i,:,:));
    cpt_refs = 0;
    cpr_domain = 0;
    for y=1:size(mask_refs,1)
        for x=1:size(mask_refs,2)
            if crt_refs(y,x)
                cpt_refs = cpt_refs + 1;
                pos_crt_refs(cpt_refs,:) = [y x];
            end
            if crt_domain(y,x)
                cpr_domain = cpr_domain + 1;
                pos_crt_domain(cpr_domain,:) = [y x];
            end
        end
    end
    
    k(i) = min(k(i), size(pos_crt_refs,1));
    crt_closest_refs = zeros(size(pos_crt_domain,1), k(i));
    for j=1:size(pos_crt_domain)
        d = sum((pos_crt_refs - pos_crt_domain(j,:)).^2,2);
        [~, idx] = mink(d, k(i));
        crt_closest_refs(j,:) = idx;
    end
    my_closest_refs{i} = crt_closest_refs;
    crt_refs = (crt_refs + squeeze(masks(i,:,:))) > 0;
end
end

