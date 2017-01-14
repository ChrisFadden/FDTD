M = csvread('./build/output/800.csv');
M = M(:,1:end-1);
h = surf(M,'EdgeColor','None','facecolor','interp'); 
colormap(parula);
colorbar;
view(2)
%exit;
