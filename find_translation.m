function tau = find_translation(R, T, kappa)
%FIND_TRANSLATION Summary of this function goes here
%   Detailed explanation goes here

[n, m] = size(R);

x0 = floor((kappa-1)*n/2)+1;
x1 = floor((kappa+1)*n/2);
y0 = floor((kappa-1)*m/2)+1;
y1 = floor((kappa+1)*m/2);

CROSS = zeros(kappa*n, kappa*m);
CROSS(x0:x1, y0:y1) = fftshift(R .* conj(T)./abs(R.*conj(T)));
cross = abs(ifft2(ifftshift(CROSS)));

ma = max(cross(:));
[ty, tx] = find(cross==ma,1,'first');

tau = [mod((ty-1)/kappa + n/2,n)-n/2 mod((tx-1)/kappa + m/2,m)-m/2];
end