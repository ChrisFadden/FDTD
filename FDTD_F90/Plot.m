fprintf('This assumes that the MaxTime is an integer multiple of 100\n');
fprintf('If you are getting missing file errors or the output looks wrong, consider changing the file read');
file = './build/output/199.csv';
M = csvread(file);
M = M(:,1:end-1);
h = surf(M,'EdgeColor','None','facecolor','interp');
colormap(jet);
colorbar;
view(2)
