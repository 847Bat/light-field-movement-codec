for /L %%F in (1 2 3 4 9 10) do (
ffmpeg -r 30 -f concat -safe 0 -i YOURPATH\04_Lists\I%%F_list.txt -s 626x434 -framerate 30 -c:v rawvideo -pix_fmt yuv422p10le YOURPATH\I%%F_13x13.yuv
)