function my_blocks = segmentation(grey_LF)
%SEGMENTATION Summary of this function goes here

resized = imresize(abs(squeeze(grey_LF(7,2,:,:)) - squeeze(grey_LF(7,12,:,:))) + ...
    abs(squeeze(grey_LF(2,7,:,:)) - squeeze(grey_LF(12,7,:,:))), [512 512]);

S = qtdecomp(resized, 0.9, [32 128]);

my_blocks = zeros(length(find(S>0)), 2, 2);
resizer = [(size(grey_LF,3)-1)/511 0 ; 0 (size(grey_LF,4)-1)/511];
counter = 1;
for i=1:512
    for j=1:512
        if S(i,j) > 0
            block_512 = [i i + S(i,j) - 1; j j + S(i,j) - 1];
            block_resized = floor(resizer*(block_512 - 1) + 1);
            if floor(resizer(1,1)*(i + S(i,j) - 2) + 1) == floor(resizer(1,1)*(i + S(i,j) - 1) + 1) && i + S(i,j) - 1 < size(grey_LF,3)
                block_resized(1,2) = block_resized(1,2) - 1;
            end
            if floor(resizer(2,2)*(j + S(i,j) - 2) + 1) < floor(resizer(2,2)*(j + S(i,j) - 1) + 1) - 1 && j + S(i,j) - 1 < size(grey_LF,4)
                block_resized(2,2) = block_resized(2,2) + 1;
            end
            my_blocks(counter, :, :) = block_resized;
            counter = counter + 1;
        end
    end
end

end

