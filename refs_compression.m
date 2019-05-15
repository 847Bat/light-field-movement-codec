function refs_10b = refs_compression(img, mask, fmt, qp)
%REFS_COMPRESSION Summary of this function goes here
%   Detailed explanation goes here

mask_4 = zeros(4, fix(numel(mask)/4 + 1));
mask_4(1:numel(mask)) = mask(:);

mask_name = '';
for i=1:size(mask_4,2)
    mask_name = [mask_name dec2hex(bin2dec(int2str(mask_4(:,i)).'))];
end
    
global_folder = ['refs/' img '/' mask_name];
mid_folder = ['refs/' img '/' mask_name '/' fmt];
folder = ['refs/' img '/' mask_name '/' fmt '/' int2str(qp)];
fname = [global_folder '/ref_list.txt'];

% Write uncompressed refs as a txt video
if ~exist(global_folder, 'dir')
    mkdir(global_folder);
    
    fileId = fopen(fname, 'w');
    
    for i=1:size(mask,1)
        for j=1:size(mask,2)
            if mask(i,j)
                im_name = sprintf('../%03d_%03d.ppm',i,j);
                fprintf(fileId, 'file ''%s'' \nduration 1\n', im_name);
            end
        end
    end

    fclose(fileId);
end

% Compress according to the parameters
if ~exist(mid_folder, 'dir')
    mkdir(mid_folder);
    setenv('PATH', getenv('PATH'));
    cmd = ['ffmpeg -r 30 -f concat -safe 0 -i ./' fname ...
        ' -s 626x434 -framerate 30 -c:v rawvideo -pix_fmt yuv' fmt ...
        'p10le ./' mid_folder '/ref_video.yuv'];
    system(cmd);
end
   
if ~exist(folder, 'dir')
    mkdir(folder);
    cmd = sprintf(['../HM-15.0+RExt-8.1/bin/TAppEncoderStatic ' ...
        '-c ../HM-15.0+RExt-8.1/cfg/encoder_lowdelay_main10.cfg ' ...
        '--InputFile=' mid_folder '/ref_video.yuv --InputBitDepth=10 ' ...
        '--Profile=main-RExt --SourceWidth=626 --SourceHeight=434 ' ...
        '--ConformanceWindowMode=1 --FrameRate=30 --FramesToBeEncoded=%d ' ...
        '--InputChromaFormat=' fmt ' --QP=%d' ' --BitstreamFile=' ...
        folder '/refs.bin --ReconFile=' folder '/refs.yuv'], ...
        sum(sum(mask)), qp);
    system(cmd);
    
    cmd = ['ffmpeg -f rawvideo -s 626x434 -r 30 -pix_fmt yuv' fmt ...
        'p10le -i ' folder '/refs.yuv ' folder '/%03d.ppm'];
    system(cmd);
end

refs = zeros(sum(sum(mask)), 434, 626, 3, 'uint16');
for i=1:sum(sum(mask))
    fname = sprintf([folder '/%03d.ppm'], i);
    refs(i,:,:,:) = imread(fname);
end

refs_10b = uint16(double(refs)/65535*1023);
end

