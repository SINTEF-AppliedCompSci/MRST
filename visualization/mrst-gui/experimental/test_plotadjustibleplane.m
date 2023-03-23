mrstModule add internal/mrst-gui spe10
layers = 1:60;
G = computeGeometry(cartGrid([60 220 numel(layers)]));
rock = SPE10_rock(layers);
%%

close all
% plotCellData(G, rock.poro, 'FaceAlpha', .3);
plotGrid(G, 'facea', 0, 'edgea', 0.05)
% plotCellData(G, rock.poro, 'FaceAlpha', .3);
axis tight off
view(30, 30)
plotAdjustiblePlane(G, log(rock.perm(:, 1)))
% fastRotateButton
%%
close all;
plotToolbar(G, rock)