fprintf('This assumes that the MaxTime is an integer multiple of 100\n');
fprintf('If you are getting missing file errors or the output looks wrong, consider changing the files read\n'); 
  set(0, 'defaultfigurevisible','off')
  
  fp = fopen('./build/output/limits.txt','r');
  limits = fscanf(fp,'%f');  
  i = 1;
 for t = 100:199   
    %Make ColorMap
    num = num2str(t);
    file = strcat('./build/output/',num,'.csv');
    M = csvread(file);
    M = M(:,1:end-1);
    h = surf(M,'EdgeColor','None','facecolor','interp'); 
    colormap(jet);
    caxis([limits(1),limits(2)]);
    colorbar;
    view(2)
 
   mov(i) = getframe(gcf);
   i = i + 1;
 end
 
 movie2avi(mov,'FDTD_Movie.avi');

exit;
