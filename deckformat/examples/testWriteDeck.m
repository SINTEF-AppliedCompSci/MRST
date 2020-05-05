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
