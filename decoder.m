function decoded_LF_10b = decoder(filename)
%DECODER Decode the encoded light field image
%   RETURN
%
% decoded_LF_10b: decoded LF picture, 4D in 10 bits and RGB

% Decompress values
disp("Unpacking data");

load([filename '_params.mat']);
load([filename '_coeffs.mat']);
load([filename '_residuals.mat']);

blocks = reshape(bi2de(bwunpack(blocks_c, prod(blocks_c_p))), blocks_c_p);
masks = reshape(bwunpack(masks_c, prod(masks_c_p)), masks_c_p);
mask_refs = reshape(bwunpack(mask_refs_c, prod(mask_refs_c_p)), mask_refs_c_p);

% Getting refs
mask_4 = zeros(4, fix(numel(mask_refs)/4 + 1));
mask_4(1:numel(mask_refs)) = mask_refs(:);

mask_name = '';
for i=1:size(mask_4,2)
    mask_name = [mask_name dec2hex(bin2dec(int2str(mask_4(:,i)).'))];
end

folder = ['refs/' img '/' mask_name '/' fmt '/' int2str(qp)];

refs_16b = zeros(sum(sum(mask_refs)), sizes_c{1}(2), sizes_c{1}(3), 3, 'uint16');

for i=1:sum(sum(mask_refs))
    fname = sprintf([folder '/%03d.ppm'], i);
    refs_16b(i,:,:,:) = imread(fname);
end
refs_10b = uint16(double(refs_16b)/65535*1023);   % Go from 16b to 10b

refs = zeros(size(refs_10b));
for i=1:size(refs,1)
    refs(i,:,:,:) = double(rgb2ycbcr10bit(squeeze(refs_10b(i,:,:,:))))/1023;
end
% Other parameters

mask_refs_tot = repmat(mask_refs, [1 1 sizes_c{1}(2) sizes_c{1}(3) 3]);
masks_tot = repmat(masks, [1 1 1 sizes_c{1}(2) sizes_c{1}(3) 3]);
decoded_LF = zeros(size(mask_refs_tot));
decoded_LF(mask_refs_tot) = refs(:);

my_closest_refs_y = closest_refs(mask_refs, masks, squeeze(nb_neighbors(:,1)));
my_closest_refs_u = closest_refs(mask_refs, masks, squeeze(nb_neighbors(:,2)));
my_closest_refs_v = closest_refs(mask_refs, masks, squeeze(nb_neighbors(:,3)));

for i=1:size(masks,1)
    disp("LEVEL " + int2str(i));
    
    taus = reshape(bi2de(bwunpack(taus_c{i}, prod(taus_c_p{i}))), taus_c_p{i}) + cast(taus_c_p2{i}, 'double');
    coeffs = reshape(bi2de(bwunpack(coeffs_c{i}, prod(coeffs_c_p{i}))), coeffs_c_p{i})*coeffs_c_p3{i}/(2^16-1) + coeffs_c_p2{i};
    
    coeffs_u = zeros(size(coeffs));
    coeffs_u(:) = 1/size(my_closest_refs_u{i},2);
    coeffs_v = zeros(size(coeffs));
    coeffs_v(:) = 1/size(my_closest_refs_v{i},2);
    
    % Predictor decoder
    disp("Predictor decoder");
    predicted_y = predictor_decoder(squeeze(refs(:,:,:,1)), taus, coeffs, blocks, my_closest_refs_y{i});
    predicted_cb = predictor_decoder(squeeze(refs(:,:,:,2)), taus, coeffs_u, blocks, my_closest_refs_u{i});
    predicted_cr = predictor_decoder(squeeze(refs(:,:,:,3)), taus, coeffs_v, blocks, my_closest_refs_v{i});
    
    % Residual decompression
    if nb_components(i) ~= 0
        disp("Residuals decompression");
        sq = reshape(bi2de(bwunpack(sq_c{i}, prod(sq_c_p2{i}))), sq_c_p2{i});
        s_sign = reshape(bwunpack(sq_c_p{i}, sq_c_p2{i}(1)), sq_c_p2{i});
        s_range = sq_c_p3{i};
        cq = reshape(bi2de(bwunpack(cq_c{i}, prod(cq_c_p2{i}))), cq_c_p2{i});
        c_sign = reshape(bwunpack(cq_c_p{i}, cq_c_p2{i}(1)), cq_c_p2{i});
        c_range = cq_c_p3{i};
        sizes = sizes_c{i};
        muq = bi2de(bwunpack(muq_c{i}, sizes(1))).' - 2^15;
        muq_range = muq_c_p{i};
        res_reconstructed = residual_decompression(sq, s_sign, s_range, cq, c_sign,...
         c_range, muq, muq_range, 4, 'log', sizes);
 
        % Getting ready for another turn
        predicted_y = max(min(predicted_y + res_reconstructed, 1),0);
    else
        predicted_y = max(min(predicted_y, 1),0);
    end
    
    predicted_cb = max(min(predicted_cb, 1),0);
    predicted_cr = max(min(predicted_cr, 1),0);
    refs_obtained = cat(4, predicted_y, predicted_cb, predicted_cr);
    refs = cat(1, refs_obtained, refs);
    
    % Stock that
    decoded_LF(squeeze(masks_tot(i,:,:,:,:,:))) = refs_obtained(:);
end

decoded_LF_10b = zeros(size(decoded_LF));
for i=1:size(decoded_LF,1)
    for j=1:size(decoded_LF,2)
        decoded_LF_10b(i,j,:,:,:) = ycbcr2rgb10bit(uint16(squeeze(decoded_LF(i,j,:,:,:))*1023));
    end
end

disp("Decoding completed.");

end

