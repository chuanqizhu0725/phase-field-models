g++ -fopenmp \
-I /opt/X11/include -L /opt/X11/lib -lX11 \
main.cpp -o main \
&& rm -f \
data/phi/*.vtk \
data/phi/*.csv \
data/con/*.csv \
data/con/*.vtk \
data/temp/*.csv \
data/interface/*.csv \
data/fraction/*.csv \
figures/phi/*.png \
figures/con/*.png \
figures/temp/*.png \
&& ./main \
&& rm main \
&& python plot2d.py