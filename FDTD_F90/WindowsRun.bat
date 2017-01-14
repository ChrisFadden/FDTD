cd build
mingw32-make
cd ..

java -cp ./gui/ GuiMain

.\build\FDTD

matlab -nosplash -nodesktop -r "run PlotMovie"

del grid.txt
