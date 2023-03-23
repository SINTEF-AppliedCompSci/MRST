%% Exercise 3.1.1
% A brick with graded grid and two cylindrical cut-outs
dx = .04;
x = 1-0.5*cos((-1:dx:1)*pi);
x = -1 - 1.5*dx+dx*cumsum(x);
y = 0:0.05:.5;
z = 0:0.04:1;
G = computeGeometry(tensorGrid(x, y, z));
cyl = @(x) x(:,1).^2 + x(:,3).^2;
G = removeCells(G,cyl(bsxfun(@minus, G.cells.centroids,[0 0 0]))<0.16);
G = removeCells(G,cyl(bsxfun(@minus, G.cells.centroids,[0 0 1]))<0.16);
clf, plotCellData(G,G.cells.volumes,'EdgeAlpha',.2); view(-20,25); axis equal tight

%% Exercise 3.1.2a
% The sum of two metaballs
g = @(x,y) (1 - min(sum(x.^2/y^2,2), 1)).^4;
f = @(x) g(bsxfun(@minus,x,[-.4,-.4, 0]),1) ...
       + g(bsxfun(@minus,x,[ .4, .4, 0]),1);
G = cartGrid([50 50 25],[2 2 1]);
G.nodes.coords = bsxfun(@minus,G.nodes.coords,[1 1 .5]);
G = computeGeometry(G);
G = removeCells(G,f(G.cells.centroids)<.3);

clf; plotGrid(G,'EdgeAlpha',.3,'FaceColor',[.6 1 1]); 
view([40,30]); axis equal; box on

%% Ecercise 3.1.2b
% The intersection of two metaballs
G = cartGrid([50 50 25],[1 1 1]);
G.nodes.coords = bsxfun(@minus,G.nodes.coords,[.5 .5 .5]);
G = computeGeometry(G);
G = removeCells(G,g(G.cells.centroids,1)<.3);
G = removeCells(G,g(bsxfun(@minus,G.cells.centroids,[-.5 0 0]),sqrt(.7))>.3);

clf; plotGrid(G,'EdgeAlpha',.3,'FaceColor',[1 .6 1]);
view([-50,30]); axis equal; box on

%% Exercise 3.1.3
% Grid form the penny data set
load penny
n = 5;
G = cartGrid([size(P),n]);
G = removeCells(G,repmat(P(:)<20,n,1));
rock.perm = repmat((950*P(:) + 50)*milli*darcy,n,1);
rock.perm = rock.perm(G.cells.indexMap);
clf, plotCellData(G,log10(rock.perm),'EdgeAlpha',.1);
view(80,80); axis tight off

%%
% A famous mouse grid
img = imread(fullfile('img','mickey.png')); img=double(img);
G = cartGrid(size(img));
img = img(:)/max(img(:));
G = removeCells(G,img==1);
rock.perm = (950*img+50)*milli*darcy; 
rock.perm=rock.perm(G.cells.indexMap);
clf, plotCellData(G,rock.perm,'EdgeAlpha',.3); view(53,80); axis tight off

%%
% Version with only n cells in vertical direction
n = 10;
img = imread(fullfile('img','mickey.png')); img=double(img(:,:,1));
G = cartGrid([size(img),n]);
img = img(:)/max(img(:));
G = removeCells(G,repmat(img==1,n,1));
rock.perm = repmat((950*img+50)*milli*darcy,n,1); 
rock.perm=rock.perm(G.cells.indexMap);
clf, plotCellData(G,rock.perm,'EdgeAlpha',.3); view(53,80); axis tight off