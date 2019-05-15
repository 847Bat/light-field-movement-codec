% Author: Irene Viola (irene.viola@epfl.ch)
% Copyright(c) Multimedia Signal Processing Group (MMSPG),
%              Ecole Polytechnique Federale de Lausanne (EPFL)
%              http://mmspg.epfl.ch
% All rights reserved.
%
% script to generate subaperture images in .ppm file format from 
% a given lenslet image (including associated lenslet image meta data)
% This script specifically adopted for the lenslet light field data set
% made available to proponents of the JPEG Pleno Lightfield Call for Proposals. 
%
%
% Usage:
%
%      saveImages( LensletImagePath, MetadataPath, OutputFolderPath )
%
% Input LensletImagePath: is the absolute or relative path of the lenslet
% image in RGB 10-bit representation that has already been demosaiced and devignetted.
%
% Input MetadataPath: is the absolute or relative path of the metadata
% corresponding to the image pointed to by LensletImagePath. The metadata
% contains the data LensletGridModel and DecodeOptions.
%
% Input OutputFolderPath: is the absolute or relative path into which the subaperture
% images will written stored (ppm files in RGB 10-bit).
%
% ---------------------------------------------------------------------------------------
%

function saveImages( LensletImagePath, MetadataPath, FolderPath )
[~, LFCol] = PPMLenslet2LF( LensletImagePath, MetadataPath );
for xx=1:15
    for yy=1:15
        
        A =squeeze(double(LFCol(xx,yy, :, :, 1:3)));
        
        A(:, end+1, :) = 0; %one-pixel padding
        
        A = A./(2^16 - 1).*(2^10 -1); %scaling to 10-bit precision
        %clipping
        A(A<0) = 0;
        A(A>1023) = 1023;
        %rounding to integer
        A = double(uint16(A));
        
        %going back to uint16
        A = A./(2^10 -1) .*(2^16-1);
        A = uint16(A);
        
        imwrite(A, sprintf('%s\\%03d_%03d.ppm', FolderPath, yy-1, xx-1), 'MaxValue', 1023);
                
    end
end