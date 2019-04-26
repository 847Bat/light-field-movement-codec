function my_psnr = lf_psnr(im,ref)
%LF_PSNR Summary of this function goes here
%   Detailed explanation goes here
[N, ~, ~] = size(im);
p = zeros(1,N);
for i=1:N
    p(i) = psnr(squeeze(im(i,:,:)),squeeze(ref(i,:,:)));

end
% disp(p);
my_psnr = mean(mean(p));
end

