function my_blocks = segmentation(grey_LF, blocks_nb)
%SEGMENTATION Summary of this function goes here

[n, m, ~, ~] = size(grey_LF);
resized = imresize(abs(squeeze(grey_LF(ceil(n/2),2,:,:)) - squeeze(grey_LF(ceil(n/2),m-1,:,:))) + ...
    abs(squeeze(grey_LF(2,ceil(m/2),:,:)) - squeeze(grey_LF(n-1,ceil(m/2),:,:))), [512 512]);

low_blocks_nb = 0.9*blocks_nb;
high_blocks_nb = 1.1*blocks_nb;
th = 0.9;
step = 0.1;
max_iter = 40;

S = qtdecomp(resized, th, [32 512]);
S2 = S;
n_iter = 0;
while length(find(S>0)) < low_blocks_nb || length(find(S>0)) > high_blocks_nb
    n_iter = n_iter + 1;
    if length(find(S>0)) < low_blocks_nb
        if length(find(S2>0)) > high_blocks_nb
            step = step / 2;
        end
        th = th - step;
    elseif length(find(S>0)) > high_blocks_nb
        if length(find(S2>0)) < low_blocks_nb
            step = step / 2;
        end
        th = th + step;
    end
    S2 = S;
    S = qtdecomp(resized, th, [32 512]);
    if n_iter > max_iter
        break
    end
    if th < 0
        th = 0.001;
    end
end

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

