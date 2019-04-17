function status = encoder(grey_LF, filename)
%ENCODER Summary of this function goes here
%   Detailed explanation goes here

% Object segmentation  (blocks)

disp("Segmentation");
resized = imresize(abs(squeeze(grey_LF(8,3,:,:)) - squeeze(grey_LF(8,13,:,:))) + ...
    abs(squeeze(grey_LF(3,8,:,:)) - squeeze(grey_LF(13,8,:,:))), [512 512]);
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
resizer = [size(grey_LF,3)/512 0 ; 0 size(grey_LF,4)/512];
counter = 1;
for i=1:512
    for j=1:512
        if S(i,j) > 0
            block_512 = [i i + S(i,j) - 1; j j + S(i,j) - 1];
            if ceil((j + S(i,j))*resizer(2,2)) - ceil((j + S(i,j) - 1)*resizer(2,2)) > 1 && j + S(i,j) < 512
                block_512(2,2) = j + S(i,j);
            end
            block_resized = ceil(resizer*block_512);
            my_blocks(counter, :, :) = block_resized;
            counter = counter + 1;
        end
    end
end

% Predictor

disp("Predictor");
% profile on
% Define reference
coded_1 = squeeze(grey_LF(2,2,:,:));
coded_2 = squeeze(grey_LF(2,6,:,:));
coded_3 = squeeze(grey_LF(2,10,:,:));
coded_4 = squeeze(grey_LF(2,14,:,:));
coded_5 = squeeze(grey_LF(6,2,:,:));
coded_6 = squeeze(grey_LF(6,6,:,:));
coded_7 = squeeze(grey_LF(6,10,:,:));
coded_8 = squeeze(grey_LF(6,14,:,:));
coded_9 = squeeze(grey_LF(10,2,:,:));
coded_10 = squeeze(grey_LF(10,6,:,:));
coded_11 = squeeze(grey_LF(10,10,:,:));
coded_12 = squeeze(grey_LF(10,14,:,:));
coded_13 = squeeze(grey_LF(14,2,:,:));
coded_14 = squeeze(grey_LF(14,6,:,:));
coded_15 = squeeze(grey_LF(14,10,:,:));
coded_16 = squeeze(grey_LF(14,14,:,:));
coded_17 = squeeze(grey_LF(4,4,:,:));
coded_18 = squeeze(grey_LF(4,8,:,:));
coded_19 = squeeze(grey_LF(4,12,:,:));
coded_20 = squeeze(grey_LF(8,4,:,:));
coded_21 = squeeze(grey_LF(8,8,:,:));
coded_22 = squeeze(grey_LF(8,12,:,:));
coded_23 = squeeze(grey_LF(12,4,:,:));
coded_24 = squeeze(grey_LF(12,8,:,:));
coded_25 = squeeze(grey_LF(12,12,:,:));
%coded_m = squeeze(grey_LF(8,8,:,:));

refs = cat(3,coded_1,coded_2, coded_3, coded_4, coded_5,coded_6, coded_7, coded_8, coded_9, coded_10, coded_11, coded_12, coded_13, coded_14, coded_15, coded_16, coded_17, coded_18, coded_19, coded_20, coded_21, coded_22, coded_23, coded_24, coded_25);

[taus, coeffs] = predictor_coder(grey_LF, refs, blocks);

% profile off
% profsave


% Decoder
disp("Decoder");
predicted = predictor_decoder(refs, taus, coeffs, blocks);

% Residual compression
disp("Residual compression");
nb_bits = 4;
nb_components = 5;
[sq, s_range, s_sign, cq, c_range, c_sign, muq, mu_range]  = ...
    residual_compression(grey_LF, predicted, nb_components, nb_bits, 'log');

% Check size of compressed values

disp("Saving data...")
% compressed_sizes(1) = numel(blocks)*10;     % size of picture < 1024
% compressed_sizes(2) = numel(refs)*16;       % refs in 2^16
% compressed_sizes(3) = numel(taus)*10;  % taus are int
% compressed_sizes(4) = numel(coeffs)*16;     % coeffs in 2^16
% compressed_sizes(5) = numel(sq)*(nb_bits+1) + 16;   % sign and max
% compressed_sizes(6) = numel(cq)*(nb_bits+1) + 16;   % sign and max
% compressed_sizes(7) = numel(muq)*16 + 16;
% compressed_size = sum(compressed_sizes);
% 
% disp("Size compressed : " + num2str(compressed_size));
% disp("Detail :");
% disp(compressed_sizes./compressed_size);
% 
% disp("bpp is : " + num2str(compressed_size/numel(grey_LF)));

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
muq_c_p = mu_range;
sizes_c = size(grey_LF);

save(filename,'blocks_c','blocks_c_p','refs_c','refs_c_p', ...
    'taus_c','taus_c_p','taus_c_p2','coeffs_c','coeffs_c_p', ...
    'coeffs_c_p2','coeffs_c_p3','sq_c','sq_c_p','sq_c_p2', ...
    'sq_c_p3', 'cq_c','cq_c_p','cq_c_p2', 'cq_c_p3', 'muq_c', ...
    'muq_c_p', 'sizes_c');

status = 0;
disp("Compression completed.");

end

