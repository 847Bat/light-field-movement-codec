for /L %%F in (1 2 3 4 9 10) do (
md YOURPATH\I%%F_HEVC
for %%B in (0,1,51) do (
x265 --input YOURPATH\I%%F_13x13.yuv --input-depth 10 --input-csp i422 --fps 30 --input-res 626x434 --output YOURPATH\I%%F_HEVC\I%%F_crf%%B.hevc --output-depth 10 --profile main422-10 --crf %%B
)
)