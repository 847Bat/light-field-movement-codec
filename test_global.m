%%
addpath(genpath('modules/subaxis'));
addpath(genpath('modules/light-field-toolbox'));
addpath(genpath('modules/light-field-graph-codec'));

load 'data/4DLF/People/Fountain_&_Vincent_2.mat'    %41.36 // 42.6 with backward (20c/5c) (41.5 without refs)
% load 'data/4DLF/Grids/Danger_de_Mort.mat'    %41.42 // 42.9 with backward (20c/5c) (41.3 without refs)

filename = 'compressed.mat';

% Transform to double
LF_double = im2double(LF).^(1/2.2);
grey_LF = 0.299 * LF_double(:,:,:,:,1) + 0.587 * LF_double(:,:,:,:,2) + 0.114 * LF_double(:,:,:,:,3);
grey_LF = floor(grey_LF*(2^16))./2^16;
grey_LF = grey_LF(2:14,2:14,:,:);


%% Encoding
status = encoder(grey_LF, filename);

%% Display efficiency
file = dir(filename);
disp("Final bpp : " + num2str(file.bytes*8/numel(grey_LF)));
stats = whos('-file',filename);
bpp = num2cell([stats.bytes].'*8/numel(grey_LF));
[stats.bpp] = bpp{:};

for i = 1:length(stats)
    if stats(i).bpp > 1e-5
        msg1 = sprintf('%s',stats(i).name);
        tab = repmat(sprintf('\t'), 1, 2 - fix(length(msg1)/8));
        msg2 =sprintf('%.5f\n',stats(i).bpp);
        fprintf([msg1, tab, msg2]);
    end
end

%% Decoding
[a,b,c,decoded_LF] = decoder(filename);

%% Assess quality
disp("mean psnr :");
disp(lf_psnr(decoded_LF, grey_LF));

%%
LFDispVidCirc(repmat(decoded_LF, [1 1 1 1 3]));

%%
disp("Unpacking data");
load(filename);

blocks = reshape(bi2de(bwunpack(blocks_c, prod(blocks_c_p))), blocks_c_p);
refs = reshape(bi2de(bwunpack(refs_c, prod(refs_c_p))), refs_c_p)/2^16;
taus = reshape(bi2de(bwunpack(taus_c, prod(taus_c_p))), taus_c_p) + cast(taus_c_p2, 'double');
coeffs = reshape(bi2de(bwunpack(coeffs_c, prod(coeffs_c_p))), coeffs_c_p)*coeffs_c_p3/(2^16-1) + coeffs_c_p2;
sq = reshape(bi2de(bwunpack(sq_c, prod(sq_c_p2))), sq_c_p2);
s_sign = reshape(bwunpack(sq_c_p, sq_c_p2(1)), sq_c_p2);
s_range = sq_c_p3;
cq = reshape(bi2de(bwunpack(cq_c, prod(cq_c_p2))), cq_c_p2);
c_sign = reshape(bwunpack(cq_c_p, cq_c_p2(1)), cq_c_p2);
c_range = cq_c_p3;
sizes = sizes_c;
muq = bi2de(bwunpack(muq_c, sizes(1)*sizes(2))).' - 2^15;
muq_range = muq_c_p;
sizes2 = sizes2_c;
taus2 = reshape(bi2de(bwunpack(taus2_c, prod(taus2_c_p))), taus2_c_p) + cast(taus2_c_p2, 'double');
coeffs2 = reshape(bi2de(bwunpack(coeffs2_c, prod(coeffs2_c_p))), coeffs2_c_p)*coeffs2_c_p3/(2^16-1) + coeffs2_c_p2;
sq2 = reshape(bi2de(bwunpack(sq2_c, prod(sq2_c_p2))), sq2_c_p2);
s_sign2 = reshape(bwunpack(sq2_c_p, sq2_c_p2(1)), sq2_c_p2);
s_range2 = sq2_c_p3;
cq2 = reshape(bi2de(bwunpack(cq2_c, prod(cq2_c_p2))), cq2_c_p2);
c_sign2 = reshape(bwunpack(cq2_c_p, cq2_c_p2(1)), cq2_c_p2);
c_range2 = cq2_c_p3;
muq2 = bi2de(bwunpack(muq2_c, sizes2(1)*sizes2(2))).' - 2^15;
muq_range2 = muq2_c_p;
%%

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
