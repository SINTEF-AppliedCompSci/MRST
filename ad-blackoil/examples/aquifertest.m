%% load modules
mrstModule add mrst-gui deckformat ad-props ad-core ad-blackoil 

mrstVerbose true
gravity on
close all

%% run two test cases

datadir = getDatasetPath('aquifertest', 'download', true);


%%

for icase = 1 : 2
    
    fnroot = sprintf('2D_OILWATER_AQ_ECLIPSE_%d', icase);
    fn = fullfile(datadir, [fnroot, '.DATA']);
    
    deck = convertDeckUnits(readEclipseDeck(fn));
    [state0, model, schedule] = initEclipseProblemAD(deck);
    
    [wellSols, states, schedulereport] = simulateScheduleAD(state0, model, schedule);
    statesEcl = loadAquiferEclipseResult(datadir, fnroot);


    %% plots
    G = model.G;
    figure
    plotToolbar(G, states);
    view([2, 2]);
    colorbar
    title(sprintf('MRST computation - Case %d', icase));
    figure
    plotToolbar(G, statesEcl);
    view([2, 2]);
    colorbar
    title(sprintf('Eclipse computation - Case %d', icase));
    
    p = states{end}.pressure;
    pecl = statesEcl{end}.pressure;
    pecl = pecl*barsa;
    err = (p - pecl)./pecl;
    figure
    plotCellData(G, err);
    view([2, 2]);
    colorbar
    title(sprintf('Relative difference - Case %d', icase));
    
end

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
