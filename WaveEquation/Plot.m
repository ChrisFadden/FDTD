  
  set(0, 'defaultfigurevisible','off');
  
  fp = fopen('output/Limits.txt','r');
  limits = fscanf(fp,'%f');  
  
  for t = 0:99   
    %Make ColorMap
    num = num2str(t);
    file = strcat('output/',num,'.csv');
    M = csvread(file);
    M = M(:,1:end-1);
    h = surf(M,'EdgeColor','None','facecolor','interp'); 
    colormap(jet);
    caxis([limits(1),limits(2)]);
    colorbar;
    view(2);
 
    mov(t+1) = getframe(gcf);
  end
  
  movie2avi(mov,'WaveMovie.avi');

exit;
