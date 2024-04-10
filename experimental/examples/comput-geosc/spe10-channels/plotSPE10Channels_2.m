plotIx = 2:4;
[wsDG, wsDGReorder, statesDG, statesDGReorder, reportsDG, reportsDGReorder] ...
                                             = deal(cell(numel(plotIx),1));
for dNo = plotIx
    [wsDG{dNo}         , statesDG{dNo}         , reportsDG{dNo}         ] ...
                          = getPackedSimulatorOutput(problems{dNo});
    [wsDGReorder{dNo}, statesDGReorder{dNo}, reportsDGReorder{dNo}] ...
                          = getPackedSimulatorOutput(reorderProblems{dNo});
end

%%

for dNo = plotIx
    st  = statesDG{dNo};
    rep = reportsDG{dNo};
    for sNo = 1:numel(st)
        st{sNo}.under = st{sNo}.s(:,1) < 0.2;
        st{sNo}.over  = st{sNo}.s(:,1) > 0.8;
    end
    statesDG{dNo} = st;
    
    st  = statesDGReorder{dNo};
    rep = reportsDGReorder{dNo};
    for sNo = 1:numel(st)
        st{sNo}.under = st{sNo}.s(:,1) < 0.2;
        st{sNo}.over  = st{sNo}.s(:,1) > 0.8;
    end
    statesDGReorder{dNo} = st;
end

%%

names = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree(plotIx)), 'unif', false);
dx = 300;
dy = 400;
pos = [-1900, 600, dx, dy];

%%

plotWellSols({wsDG{:}, wsDGReorder{:}})
legend(names);

%%
close all

for dNo = plotIx
    if ~isempty(statesDG{dNo})
        figure('name', names{dNo}, 'Position', pos+[dx*(dNo-1),0,0,0])
        plotToolbar(G, statesDG{dNo})
        axis equal tight
        colormap(jet)
    end
    if ~isempty(statesDGReorder{dNo})
        figure('name', [names{dNo}, ' reorder'], 'Position', pos+[dx*(dNo-1),-dy*1.4,0,0])
        plotToolbar(G, statesDGReorder{dNo})
        axis equal tight
        colormap(jet)
    end
end

%%

close all
sNo = numel(statesDG{end});
for dNo = plotIx
    figure('name', names{dNo}, 'Position', pos+[dx*(dNo-1),0,0,0])
    plotCellData(G, statesDG{dNo}{sNo}.s(:,2), 'edgec', 'none');
    caxis([0.2,0.8]);
    axis equal tight
    colormap(jet)
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
