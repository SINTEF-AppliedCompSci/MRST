clc; clear all; close all;

load('showInnerOuterCellsElephant.mat')
plotGrid(Gc,'faceAlpha', .2)
axis equal off

w = 0;
h = 2;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
print(gcf, '-dpdf', '../../tex/thesis/fig/Elephant.pdf');
