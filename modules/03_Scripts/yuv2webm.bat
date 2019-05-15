for /L %%F in (1 2 3 4 9 10) do (
md YOURPATH\I%%F_vp9
for %%B in (10788 718 144 20) do (
vpxenc --i422 --input-bit-depth=10 --profile=3 -w 626 -h 434 --target-bitrate=%%B --cq-level=0 --bit-depth=10 --codec=vp9 --fps=30000/1000 --best -o I%%F_vp9/I%%F_%%B.webm I%%F_13x13.yuv
)
)