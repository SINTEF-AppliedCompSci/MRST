function state = Initialize(model)
    Swi        = model.experiment.rock.Swi.value;
    pini       = model.experiment.schedule.pini.value;
    simulation = model.simulation;
    G          = model.grid.G;
    state      = initState(G, [], pini, [Swi 1-Swi]); % filled with phase 2 (oil)
    if(isfield(simulation,'bCells'))
        firstCellSw    = simulation.bCells.firstCellSw;
        lastCellSw     = simulation.bCells.lastCellSw;
        state.s(1,1)   = firstCellSw.value;
        state.s(1,2)   = 1 - firstCellSw.value;           
        state.s(end,1) = lastCellSw.value;
        state.s(end,2) = 1 - lastCellSw.value;
    end   
end