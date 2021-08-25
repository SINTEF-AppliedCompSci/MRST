deck = convertDeckUnits(readEclipseDeck('HISTSW.DATA'));

G = computeGeometry(initEclipseGrid(deck));

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);


model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false);

model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true);

%model.FacilityModel = UniformFacilityModel(model);


%%
W = processWells(G, rock, deck.SCHEDULE.control(end)); % Pick arbitrary control which is not empty

for i = 1:numel(W)
    if ~strcmp(W(i).type,'rate')
        W(i).type = 'bhp';
        W(i).val = 100*barsa;
        W(i).compi = [0 0];
        W(i).status = true;
        W(i).cstatus(W(i).cstatus==0) = true;
        W(i).lims = [];
    else
        W(i).val = 0.1;
        W(i).compi =  [1 0];
        W(i).status = true;
        W(i).cstatus(W(i).cstatus==0) = true;
        W(i).lims = [];
    end
end

dT = rampupTimesteps(10*year, 30*day);
schedule = simpleSchedule(dT, 'W', W);
%%
state0 = initStateDeck(model, deck);

nonlinear = NonLinearSolver();
nonlinear.LinearSolver = selectLinearSolverAD(model);

%%
 [wellSols, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nonlinear);
%% 
 plotWellSols(wellSols, cumsum(dT))
 
 
 