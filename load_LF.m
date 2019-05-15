function LF = load_LF(img)
%LOAD_LF Summary of this function goes here
%   Detailed explanation goes here

folder = ['refs/' img];
LF = zeros(15, 15, 434, 626, 3, 'uint16');

for i=0:14
    for j=0:14
        imname = sprintf('%03d_%03d.ppm', j, i);
        LF(i+1,j+1,:,:,:) = imread([folder '/' imname]);
    end
end

end

