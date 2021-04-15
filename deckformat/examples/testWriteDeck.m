%% Write ELCIPSE/OPM-type deck to files
% This example demonstrates the use of model2deck/writeDeck to create/wtite 
% ECLIPSE/OPM-type input-decks from deck-structures and/or MRST-models

mrstModule add ad-core

mrstVerbose true 

%% Read and write deck for SPE1
pth = getDatasetPath('spe1');
fn  = fullfile(pth, 'BENCH_SPE1.DATA');
deck_orig = convertDeckUnits( readEclipseDeck(fn));

% write a new deck using metric units and list newly created files
outpth = fullfile(mrstDataDirectory(), 'SPE1_METRIC') ;
writeDeck(deck_orig, outpth, 'unit', 'metric');
fprintf('\nContents of %s:\n', outpth)
dir(outpth)

%% Read new deck and print contents of DATA-file
fn_new = fullfile(outpth, 'SPE1_METRIC.DATA');
deck_new = convertDeckUnits( readEclipseDeck(fn_new));
% view contents of newly created of DATA-file
mstat = get(0, 'more');
more('on');
type(fn_new); 
more(mstat);

%% Check that the two decks indeed produce the same results
[initState, model, schedule] = initEclipseProblemAD(deck_orig);
ws1 = simulateScheduleAD(initState, model, schedule);
[initState2, model2, schedule2] = initEclipseProblemAD(deck_new);
ws2 = simulateScheduleAD(initState2, model2, schedule2);

plotWellSols({ws1, ws2}, {cumsum(schedule.step.val), cumsum(schedule2.step.val)}, ...
             'datasetnames', {'Original deck', 'Generated deck'});
    
%% Create deck from MRST-model        
% Create a second deck directly from the MRST model. The grid is represented 
% as 1D with all connections given as NNCs and all transmissibilities of
% the 1D-grid set to zero
%
% The function model2deck does not support converting fluid-structures to
% PROPS, so this is copied from the original deck

deck_mrst = model2Deck(model, schedule, 'deck', deck_orig);
% view contents of GRID-section
fprintf('\nContents of GRID-section:\n')
disp(deck_mrst.GRID)
% view contents of EDIT-section
fprintf('Contents of EDIT-section:\n')
disp(deck_mrst.EDIT)

% As above, deck_mrst can written to files by writeDeck and run with other
% appropriate simulators. MRST does currently not support initializing a
% model based on such a deck due to lack of required geometry information
% (likely to change soon)

%% Copyright Notice
%
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
