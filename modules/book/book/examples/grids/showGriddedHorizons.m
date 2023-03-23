%%
[y,x,z]  = peaks(15); z = z+8;
horizons = {struct('x',x,'y',y,'z',z),struct('x',x,'y',y,'z',z+8)};
grdecl   = convertHorizonsToGrid(horizons,'layers', 4);
G        = processGRDECL(grdecl);
figure, plotGrid(G); view(3); axis tight  off

figure
surf(horizons{1}.x,horizons{1}.y,horizons{1}.z, 'EdgeC','r','FaceC',[.8 .8 .8]),  hold on
mesh(horizons{2}.x,horizons{2}.y,horizons{2}.z, 'EdgeC','b','FaceC',[.7 .7 .7])
plot3([horizons{1}.x(:) horizons{2}.x(:)]',...
    [horizons{1}.y(:) horizons{2}.y(:)]',...
    [horizons{1}.z(:) horizons{2}.z(:)]',':k');
axis tight off, set(gca,'ZDir','reverse');
set(gca,'XLim',[-3.05 3.05]);


%%
[n,m] = deal(30);
horizons = {struct('x',x,'y',y,'z',z),struct('x',x+.5,'y',2*y+1,'z',z+10)};
grdecl = convertHorizonsToGrid(horizons,'dims',[n m], 'layers', 3);
G = processGRDECL(grdecl);

figure
h1 = surf(horizons{1}.x,horizons{1}.y,horizons{1}.z-.1, ...
    'EdgeC','r','FaceC',[.8 .8 .8]);  hold on
h2 = mesh(horizons{2}.x,horizons{2}.y,horizons{2}.z, ...
    'EdgeC','b','FaceC',[.7 .7 .7]);

xmin = min(cellfun(@(h) min(h.x(:)), horizons));
xmax = max(cellfun(@(h) max(h.x(:)), horizons));
ymin = min(cellfun(@(h) min(h.y(:)), horizons));
ymax = max(cellfun(@(h) max(h.y(:)), horizons));

[xi,yi] = ndgrid(linspace(xmin,xmax,n+1), linspace(ymin,ymax,m+1));
hi = mesh(xi,yi,26*ones(size(xi)),'FaceC','none');
hg = plotGrid(G);
hb1 = plot3(...
    [-3  3  3 -3 -3 NaN -3 -3 NaN  3  3 NaN  3  3 NaN -3 -3], ...
    [-3 -3  3  3 -3 NaN -3 -3 NaN -3 -3 NaN  3  3 NaN  3  3],...
    [26 26 26 26 26 NaN 26  8 NaN 26  8 NaN 26  8 NaN 26  8],'r-','LineWidth',2);
hb2 = plot3(...
    [-2  4  4 -2 -2 NaN -2 -2 NaN  4  4 NaN  4  4 NaN -2 -2]-.5, ...
    [-5 -5  7  7 -5 NaN -5 -5 NaN -5 -5 NaN  7  7 NaN  7  7],...
    [26 26 26 26 26 NaN 26 18 NaN 26 18 NaN 26 18 NaN 26 18],'b-','LineWidth',2);
axis tight off, view(-140,30)