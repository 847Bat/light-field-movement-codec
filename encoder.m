function status = encoder(grey_LF, filename)
%ENCODER Summary of this function goes here
%   Detailed explanation goes here

% Object segmentation  (blocks)
disp("Segmentation");
blocks = segmentation(grey_LF);

blocks_c = bwpack(de2bi(blocks));
blocks_c_p = cast(size(blocks), 'uint32');

% Initializing multilevel predictor
mask_refs = zeros(size(grey_LF), 'logical');
mask_refs(1:4:end,1:4:end, :, :) = 1;
refs = reshape(grey_LF(mask_refs), [sum(sum(mask_refs(:,:,1,1))) size(grey_LF,3) size(grey_LF, 4)]);
refs_c = bwpack(de2bi(refs*(2^16)));
refs_c_p = cast(size(refs), 'uint32');
mask_refs_c = bwpack(mask_refs(:));
mask_refs_c_p = cast(size(mask_refs), 'uint32');


nb_components = [15 5 1];
    

% Building masks
masks = zeros(1, size(grey_LF,1), size(grey_LF,2), size(grey_LF,3), size(grey_LF,4), 'logical');
% masks(1, 3:4:end, 3:4:end, :, :) = 1;
% masks(2, 1:2:end, 1:2:end, :, :) = 1;
% masks(2, masks(1,:,:,:,:)) = 0;
% masks(2, mask_refs) = 0;
% masks(3, :, :, :, :) = 1;
% masks(3, masks(1,:,:,:,:)) = 0;
% masks(3, masks(2,:,:,:,:)) = 0;
% masks(3, mask_refs) = 0;
masks(1, :, :, :, :) = 1;
% masks(1, mask_refs) = 0;
% masks(2, :, :, :, :) = 1;
% masks(2, mask_refs) = 0;
masks_c = bwpack(masks(:));
masks_c_p = cast(size(masks), 'uint32');

for i=1:size(masks,1)
    disp("LEVEL " + int2str(i));
    % Predictor
    disp("Predictor");
    crt_domain = reshape(grey_LF(masks(i,:,:,:,:)), [sum(sum(masks(i,:,:,1,1))) size(grey_LF,3) size(grey_LF, 4)]);
    [taus, coeffs] = predictor_coder(crt_domain, refs, blocks);
    coeffs = floor((coeffs - min(coeffs(:)))*(2^16-1)/(max(coeffs(:))-min(coeffs(:))))*(max(coeffs(:))-min(coeffs(:)))/(2^16-1) + min(coeffs(:));

    % Decoder
    disp("Decoder");
    predicted = predictor_decoder(refs, taus, coeffs, blocks);
    save('tmp2.mat','predicted');
    disp(lf_psnr(predicted, crt_domain));

    % Residual compression
    disp("Residual compression");
    nb_bits = 4;
    [sq, s_range, s_sign, cq, c_range, c_sign, muq, muq_range]  = ...
        residual_compression(crt_domain, predicted, nb_components(i), nb_bits, 'log');
    
    % Saving data    
    taus_c{i} = bwpack(de2bi(taus(:) - min(taus(:))));
    taus_c_p{i} = cast(size(taus), 'uint32');
    taus_c_p2{i} = min(taus(:));
    coeffs_c{i} = bwpack(de2bi(floor((coeffs - min(coeffs(:)))*(2^16-1)/(max(coeffs(:))-min(coeffs(:))))));
    coeffs_c_p{i} = cast(size(coeffs), 'uint32');
    coeffs_c_p2{i} = min(coeffs(:));
    coeffs_c_p3{i} = max(coeffs(:))-min(coeffs(:));
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
    sizes_c{i} = size(crt_domain);

    % Getting ready for next level
    res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign,...
     c_range, muq, muq_range, nb_bits, 'log', size(crt_domain));

    disp(lf_psnr(predicted + res_reconstructed, crt_domain));
    
    refs_obtained = max(min(predicted + res_reconstructed, 1),0);
    refs = cat(1, refs_obtained, refs);
end


% Check size of compressed values

disp("Saving data...")

save(filename,'blocks_c','blocks_c_p','refs_c','refs_c_p', ...
    'taus_c','taus_c_p','taus_c_p2','coeffs_c','coeffs_c_p', ...
    'coeffs_c_p2','coeffs_c_p3','sq_c','sq_c_p','sq_c_p2', ...
    'sq_c_p3', 'cq_c','cq_c_p','cq_c_p2', 'cq_c_p3', 'muq_c', ...
    'muq_c_p', 'sizes_c', 'masks_c', 'masks_c_p', ...
    'mask_refs_c', 'mask_refs_c_p');

status = 0;
disp("Compression completed.");

end

