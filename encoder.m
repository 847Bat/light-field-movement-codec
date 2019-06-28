function status = encoder(LF, filename, mask_refs, masks, nb_components, qp, fmt, img, nb_neighbors, nb_blocks)
%ENCODER Encode the 4D light field image
%   
% INPUTS:
%
% LF:           4D light field encoded in 10 bits in RGB of size n1xn2 x 
%               n3xn4 (n1xn2 are the views, n3xn4 the pixels).
% filename:     Prefix of the files containing the compressed data
% mask_refs:    Logical array of size n1xn2 with 1 where the view is a
%               reference.
% masks:        Logical array of size mxn1xn2 where m is the number of
%               levels. For each level, masks(i,:,:) is the equivalent of
%               mask_refs but the 1 represent the vies to predict in this
%               level.
% nb_components:Array of size m containing the number of components for the
%               PCA encoding of the residuals at each level
% qp:           Quality parameter for the HEVC encoding
% fmt:          Format of pixel in the HEVC encoding (420, 422 or 444)
% img:          Name of the image ('I09' for instance). The files
%               %03d_%03d.ppm, encoded in 10 bits should be in the folder
%               ./refs/img
% nb_neighbors: Array of size mx3 representing the number of reference to
%               consider for each prediction, for each channel (YUV).
%               Careful, the nb for U and V has to be lower than for Y.
% nb_blocks:    Number of blocks wanted in the segmentation (final number
%               will be between 0.9*nb_blocks and 1.1*nb_blocks
%
% RETURN:
%
% status
%
% The references have been encoded in ./refs/img/fmt/qp/refs.bin (decoded
% refs are in ./refs/img/fmt/qp/ in the form %03d.ppm).
%
% The coefficients and the translations for the predictor have been encoded
% in filename_coeffs.mat
% The residuals have been encoded in filename_residuals.mat
% The remaining parameters have been encoded in filename.params.mat

% Transform LF from RGB, 10 bits to YUV, double
prec = 1;

LF_ycbcr = zeros(size(LF), 'uint16');
for i=1:size(LF,1)
    for j=1:size(LF,2)
        LF_ycbcr(i,j,:,:,:) = rgb2ycbcr10bit(squeeze(LF(i,j,:,:,:)));
    end
end
LF_double = double(LF_ycbcr)/1023;
grey_LF = squeeze(LF_double(:,:,:,:,1));

% Object segmentation  (blocks)
disp("Segmentation");
blocks = segmentation(grey_LF/max(grey_LF(:)), nb_blocks);

% Initializing multilevel predictor
% Compressing the references
refs_10b = refs_compression(img, mask_refs, fmt, qp);

% Transform the refs from RGB, 10 bits to YUV, double
refs = zeros(size(refs_10b));
for i=1:size(refs_10b,1)
    refs(i,:,:,:) = double(rgb2ycbcr10bit(squeeze(refs_10b(i,:,:,:))))/1023;
end

% Selecting only the Y channel
refs = squeeze(refs(:,:,:,1));

% Building masks
masks_grey = repmat(masks, [1 1 1 size(LF,3) size(LF,4)]);

% Getting the k closest refs for each view to predict
my_closest_refs_grey = closest_refs(mask_refs, masks, squeeze(nb_neighbors(:,1)));

% Multilevel predictor, iterate on the levels
for i=1:size(masks,1)
    disp("LEVEL " + int2str(i));
    % Predictor
    disp("Predictor");
    crt_domain = reshape(grey_LF(squeeze(masks_grey(i,:,:,:,:))), [sum(sum(masks(i,:,:))) size(grey_LF,3) size(grey_LF, 4)]);
    [taus, coeffs] = predictor_coder(crt_domain, refs, blocks, my_closest_refs_grey{i}, prec);
    coeffs = floor((coeffs - min(coeffs(:)))*(2^16-1)/(max(coeffs(:))-min(coeffs(:))))*(max(coeffs(:))-min(coeffs(:)))/(2^16-1) + min(coeffs(:));

    % Decoder
    disp("Decoder");
    predicted = predictor_decoder(refs, taus, coeffs, blocks, my_closest_refs_grey{i});
    disp(lf_psnr(predicted, crt_domain));

    % Residual compression
    disp("Residual compression");
    nb_bits = 4;
    if nb_components(i) ~= 0
        [sq, s_range, s_sign, cq, c_range, c_sign, muq, muq_range]  = ...
            residual_compression(crt_domain, predicted, nb_components(i), nb_bits, 'log');
        
        sq_c{i} = bwpack(de2bi(sq));
        sq_c_p{i} = bwpack(s_sign);
        sq_c_p2{i} = cast(size(sq), 'uint32');
        sq_c_p3{i} = s_range;
        cq_c{i} = bwpack(de2bi(cq));
        cq_c_p{i} = bwpack(c_sign);
        cq_c_p2{i} = cast(size(cq), 'uint32');
        cq_c_p3{i} = c_range;
        muq_c{i} = bwpack(de2bi(muq + 2^15));
        muq_c_p{i} = muq_range;
        
        % Getting ready for next level
        res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign,...
         c_range, muq, muq_range, nb_bits, 'log', size(crt_domain));

        refs_obtained = max(min(predicted + res_reconstructed, 1),0);
    else
        sq_c{i} = 0;
        sq_c_p{i} = 0;
        sq_c_p2{i} = 0;
        sq_c_p3{i} = 0;
        cq_c{i} = 0;
        cq_c_p{i} = 0;
        cq_c_p2{i} = 0;
        cq_c_p3{i} = 0;
        muq_c{i} = 0;
        muq_c_p{i} = 0;
        refs_obtained = max(min(predicted, 1),0);
    end
    
    % Saving data    
    save('taus.mat', 'taus');
    taus_c{i} = bwpack(de2bi((taus(:) - min(taus(:)))*prec));
    taus_c_p{i} = cast(size(taus), 'uint32');
    taus_c_p2{i} = min(taus(:));
    coeffs_c{i} = bwpack(de2bi(floor((coeffs - min(coeffs(:)))*(2^16-1)/(max(coeffs(:))-min(coeffs(:))))));
    coeffs_c_p{i} = cast(size(coeffs), 'uint32');
    coeffs_c_p2{i} = min(coeffs(:));
    coeffs_c_p3{i} = max(coeffs(:))-min(coeffs(:));
    sizes_c{i} = size(crt_domain);

    disp(lf_psnr(refs_obtained, crt_domain));
    refs = cat(1, refs_obtained, refs);
end


% Check size of compressed values

disp("Saving data...")


blocks_c = bwpack(de2bi(blocks));
blocks_c_p = cast(size(blocks), 'uint32');
save([filename '_coeffs.mat'],'blocks_c','blocks_c_p', ...
    'taus_c','taus_c_p','taus_c_p2','coeffs_c','coeffs_c_p', ...
    'coeffs_c_p2','coeffs_c_p3');

mask_refs_c = bwpack(mask_refs(:));
mask_refs_c_p = cast(size(mask_refs), 'uint32');
masks_c = bwpack(masks(:));
masks_c_p = cast(size(masks), 'uint32');
save([filename '_params.mat'], 'nb_components', 'sizes_c', 'masks_c', 'masks_c_p', ...
    'mask_refs_c', 'mask_refs_c_p', 'qp', 'img', 'fmt', 'nb_neighbors', 'prec');

save([filename '_residuals.mat'], 'sq_c','sq_c_p','sq_c_p2', ...
    'sq_c_p3', 'cq_c','cq_c_p','cq_c_p2', 'cq_c_p3', 'muq_c', ...
    'muq_c_p');

status = 0;
disp("Compression completed.");

end

