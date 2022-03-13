[wsDG0 , stDG0 ] = getPackedSimulatorOutput(problemDG0);
[wsDG1 , stDG1 ] = getPackedSimulatorOutput(problemDG1);
[wsWENO, stWENO] = getPackedSimulatorOutput(problemWENO);

%%

close all

pba  = [1,3,0.4];
azel = [45,25];
pos  = [0,0,900,600];
lpos = [1500, 9500, 0];
% lpos = [1,1,2];

cmap = winter();
cmap = cmap(end:-1:1,:);
sMin = 0.1;

figDir = fullfile(mrstPath('dg'), 'examples', 'siamgs-19', 'fig', '3ph-pebi');
savepng = @(name) print(fullfile(figDir, name), '-dpng', '-r400');
% savepng = @(name) [];

W = schedule.control(1).W;
[W.name] = deal('Inj', 'Prod');

cmap = winter();

pw = @() plotWell(G, W, 'color', 'k', 'height', 20);

%%

close all

sMin = 0.2;
steps = [30, 45, 53];

for sNo = steps
    
    if problemDG0.OutputHandlers.states.numelData >= sNo
        figure('name', 'dG(0)', 'pos', pos)
        s = stDG0{sNo}.s(:,3);
        plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
        plotCellData(G, s, s>sMin);
        pw();
        axis equal off
        caxis([sMin,0.8])
        view(azel)
        colormap(cmap)
        savepng(['3ph-pebi-sat-dg0-', num2str(sNo)]);
    end

    if problemDG1.OutputHandlers.states.numelData >= sNo
        figure('name', 'dG(1)', 'pos', pos)
        s = stDG1{sNo}.s(:,3);
        plotGrid(G, 'facec', 'none', 'edgealpha', 0.2);
        plotCellData(G, s, s>sMin);
        pw();
        axis equal off
        caxis([sMin,0.8])
        view(azel)
        colormap(cmap)
        savepng(['3ph-pebi-sat-dg1-', num2str(sNo)]);
    end
    
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
