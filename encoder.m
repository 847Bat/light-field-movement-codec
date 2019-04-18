function status = encoder(grey_LF, filename)
%ENCODER Summary of this function goes here
%   Detailed explanation goes here

% Object segmentation  (blocks)

disp("Segmentation");
resized = imresize(abs(squeeze(grey_LF(7,2,:,:)) - squeeze(grey_LF(7,12,:,:))) + ...
    abs(squeeze(grey_LF(2,7,:,:)) - squeeze(grey_LF(12,7,:,:))), [512 512]);
% resized = imresize(squeeze(grey_LF(8,8,:,:)), [512 512]);
S = qtdecomp(resized, 0.9, [32 128]);
blocks = repmat(uint8(0),size(S));

for dim = [512 256 128 64 32 16 8 4 2 1]    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;
% imshow(blocks,[])

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
blocks = my_blocks;

% Predictor

disp("Predictor");

% Define reference
refs = reshape(grey_LF(1:4:end,1:4:end, :, :), [4*4 size(grey_LF,3) size(grey_LF, 4)]);
domain1 = grey_LF(1:2:end, 1:2:end, :, :);
[taus, coeffs] = predictor_coder(domain1, refs, blocks);

% Decoder
disp("Decoder");
predicted = predictor_decoder(refs, taus, coeffs, blocks);

disp(lf_psnr(predicted, domain1));

% Residual compression
disp("Residual compression");
nb_bits = 4;
nb_components = 10;
[sq, s_range, s_sign, cq, c_range, c_sign, muq, muq_range]  = ...
    residual_compression(domain1, predicted, nb_components, nb_bits, 'log');

res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign,...
 c_range, muq, muq_range, nb_bits, 'log', size(domain1));

disp(lf_psnr(predicted + res_reconstructed, domain1));

refs2 = reshape(max(min(predicted + res_reconstructed, 1),0), ...
    [size(predicted,1)*size(predicted,2) size(predicted,3) size(predicted, 4)]);

domain2 = grey_LF;
[taus2, coeffs2] = predictor_coder(domain2, refs2, blocks);

% Decoder
disp("Decoder");
predicted2 = predictor_decoder(refs2, taus2, coeffs2, blocks);

disp(lf_psnr(predicted2, domain2));

% Residual compression
disp("Residual compression");
nb_bits = 4;
nb_components = 5;
[sq2, s_range2, s_sign2, cq2, c_range2, c_sign2, muq2, muq_range2]  = ...
    residual_compression(domain2, predicted2, nb_components, nb_bits, 'log');

% Check size of compressed values

disp("Saving data...")

blocks_c = bwpack(de2bi(blocks));
blocks_c_p = cast(size(blocks), 'uint32');
refs_c = bwpack(de2bi(refs*(2^16)));
refs_c_p = cast(size(refs), 'uint32');
taus_c = bwpack(de2bi(taus(:) - min(taus(:))));
taus_c_p = cast(size(taus), 'uint32');
taus_c_p2 = min(taus(:));
coeffs_c = bwpack(de2bi(floor((coeffs - min(coeffs(:)))*(2^16-1)/(max(coeffs(:))-min(coeffs(:))))));
coeffs_c_p = cast(size(coeffs), 'uint32');
coeffs_c_p2 = min(coeffs(:));
coeffs_c_p3 = max(coeffs(:))-min(coeffs(:));
sq_c = bwpack(de2bi(sq));
sq_c_p = bwpack(s_sign);
sq_c_p2 = cast(size(sq), 'uint32');
sq_c_p3 = s_range;
cq_c = bwpack(de2bi(cq));
cq_c_p = bwpack(c_sign);
cq_c_p2 = cast(size(cq), 'uint32');
cq_c_p3 = c_range;
muq_c = bwpack(de2bi(muq + 2^15));
muq_c_p = muq_range;
sizes_c = size(domain1);

taus2_c = bwpack(de2bi(taus2(:) - min(taus2(:))));
taus2_c_p = cast(size(taus2), 'uint32');
taus2_c_p2 = min(taus2(:));
coeffs2_c = bwpack(de2bi(floor((coeffs2 - min(coeffs2(:)))*(2^16-1)/(max(coeffs2(:))-min(coeffs2(:))))));
coeffs2_c_p = cast(size(coeffs2), 'uint32');
coeffs2_c_p2 = min(coeffs2(:));
coeffs2_c_p3 = max(coeffs2(:))-min(coeffs2(:));
sq2_c = bwpack(de2bi(sq2));
sq2_c_p = bwpack(s_sign2);
sq2_c_p2 = cast(size(sq2), 'uint32');
sq2_c_p3 = s_range2;
cq2_c = bwpack(de2bi(cq2));
cq2_c_p = bwpack(c_sign2);
cq2_c_p2 = cast(size(cq2), 'uint32');
cq2_c_p3 = c_range2;
muq2_c = bwpack(de2bi(muq2 + 2^15));
muq2_c_p = muq_range2;
sizes2_c = size(domain2);

save(filename,'blocks_c','blocks_c_p','refs_c','refs_c_p', ...
    'taus_c','taus_c_p','taus_c_p2','coeffs_c','coeffs_c_p', ...
    'coeffs_c_p2','coeffs_c_p3','sq_c','sq_c_p','sq_c_p2', ...
    'sq_c_p3', 'cq_c','cq_c_p','cq_c_p2', 'cq_c_p3', 'muq_c', ...
    'muq_c_p', 'sizes_c', ...
    'taus2_c','taus2_c_p','taus2_c_p2','coeffs2_c','coeffs2_c_p', ...
    'coeffs2_c_p2','coeffs2_c_p3','sq2_c','sq2_c_p','sq2_c_p2', ...
    'sq2_c_p3', 'cq2_c','cq2_c_p','cq2_c_p2', 'cq2_c_p3', 'muq2_c', ...
    'muq2_c_p', 'sizes2_c');

status = 0;
disp("Compression completed.");

end

