gcc main.c -o main && rm -f data/phi/*.csv figures/phi/*.png data/con/*.csv figures/con/*.png data/temp/*.csv figures/temp/*.png && ./main && python plot.py && rm main