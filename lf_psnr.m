function my_psnr = lf_psnr(im,ref)
%LF_PSNR Summary of this function goes here
%   Detailed explanation goes here
[M, N, ~, ~] = size(im);
p = zeros(M,N);
for i=1:M
    for j=1:N
        p(i,j) = psnr(squeeze(im(i,j,:,:)),squeeze(ref(i,j,:,:)));
    end
end
disp(p);
my_psnr = mean(mean(p));
end

