

plotIx = 1:3;
[wsDG, statesDG, reports] = deal(cell(numel(plotIx),1));
for dNo = plotIx
    [wsDG{dNo}, statesDG{dNo}, reports{dNo}] = getPackedSimulatorOutput(reorderProblems{dNo});
end

%%

plotWellSols(wsDG)

%%

names = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree), 'unif', false);
dx = 350;
pos = [-1900, 275, dx, 600];

%%
close all

for dNo = plotIx
    figure('name', names{dNo}, 'Position', pos+[dx*(dNo-1),0,0,0])
    plotToolbar(G, statesDG{dNo})
    axis equal tight
end

%%

close all
sNo = numel(statesDG{end});
for dNo = plotIx
    figure('name', names{dNo}, 'Position', pos+[dx*(dNo-1),0,0,0])
    plotCellData(G, statesDG{dNo}{sNo}.s(:,1), 'edgec', 'none');
    caxis([0.2,0.8]);
    axis equal tight
end
