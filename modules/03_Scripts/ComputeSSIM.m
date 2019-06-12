% Author: Irene Viola (irene.viola@epfl.ch)
% Copyright(c) Multimedia Signal Processing Group (MMSPG),
%              Ecole Polytechnique Federale de Lausanne (EPFL)
%              http://mmspg.epfl.ch
% All rights reserved.
%
% script which computes SSIM values as defined in JPEG Pleno Call for 
% Proposals on Light Field Coding.
%
% Input:
% R - LF data structure created from uncompressed lenslet image in uint16, 10-bit representation [0-1023]
% I - LF data structure created from compressed bitstream in uint16, 10-bit representation [0-1023]
%
% Please note: The correct function of this script assumes a 10-bit representation of both R and I
%
%
% Output:
%  SSIM values per channel and for YUV image according to the CfP mentioned above.
%
% Note:
% Removing the comments at the end allows to save SSIM values to SSIMs.mat in the current directory
%
% ---------------------------------------------------------------------------------------
%

function [SSIM_Y, SSIM_U, SSIM_V, SSIM_YUV, SSIM_Y_mean, SSIM_U_mean, SSIM_V_mean, SSIM_YUV_mean] = ComputeSSIM(I, R)
% removing weighting channel if needed
I = I(:,:,:,:,1:3);
R = R(:,:,:,:,1:3);
m = size(I,1);
n = size(I,2);
for k = 1:m
     for l = 1:n
        Iyuv = rgb2ycbcr10bit(squeeze(I(k,l,:,:,:)));
        Ryuv = rgb2ycbcr10bit(squeeze(R(k,l,:,:,:)));
        Iyuv = double(Iyuv)./1023;
        Ryuv = double(Ryuv)./1023;
    
        SSIM_Y(k,l) = ssim(Iyuv(:,:,1), Ryuv(:,:,1));
        SSIM_U(k,l) = ssim(Iyuv(:,:,2), Ryuv(:,:,2));
        SSIM_V(k,l) = ssim(Iyuv(:,:,3), Ryuv(:,:,3));
        
        
        SSIM_YUV(k,l) = (6*SSIM_Y(k,l)+SSIM_U(k,l)+SSIM_V(k,l))/8;
    end
end

SSIM_Y(isinf(SSIM_Y)) = NaN;
SSIM_U(isinf(SSIM_U)) = NaN;
SSIM_V(isinf(SSIM_V)) = NaN;
SSIM_YUV(isinf(SSIM_YUV)) = NaN;

SSIM_Y_mean = nanmean(SSIM_Y(:));
SSIM_U_mean = nanmean(SSIM_U(:));
SSIM_V_mean = nanmean(SSIM_V(:));
SSIM_YUV_mean = nanmean(SSIM_YUV(:));

% save('SSIMs.mat','SSIM_Y', 'SSIM_U', 'SSIM_V', 'SSIM_YUV', 'SSIM_Y_mean', ' SSIM_U_mean', 'SSIM_V_mean','SSIM_YUV_mean')
