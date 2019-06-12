function LF = load_LF(img, n, m, o, p)
%LOAD_LF Summary of this function goes here
%   Detailed explanation goes here

folder = ['refs/' img];
LF = zeros(n, m, o, p, 3, 'uint16');

for i=0:n-1
    for j=0:m-1
        imname = sprintf('%03d_%03d.ppm', j, i);
        LF(i+1,j+1,:,:,:) = imread([folder '/' imname]);
    end
end

end

