g++ -fopenmp \
-I /opt/X11/include -L /opt/X11/lib -lX11 \
main.cpp -o main \
&& rm -f \
data/con/*.vtk \
data/con/*.csv \
data/interface/*.csv \
data/fraction/*.csv \
figures/con_xy/*.png \
figures/con_yz/*.png \
&& ./main \
&& rm main