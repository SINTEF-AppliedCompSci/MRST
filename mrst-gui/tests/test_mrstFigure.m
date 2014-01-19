mrstModule add mrst-gui
G = cartGrid([10 10 1]);
%% Plot something
close all
h = mrstFigure();
plotGrid(G)

%% Plot something in background
mrstFigure(h, 'background')
hold on
plot([0 10], [0 10])