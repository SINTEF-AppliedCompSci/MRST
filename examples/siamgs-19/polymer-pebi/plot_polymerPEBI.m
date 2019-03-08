[wsDG0, stDG0] = getPackedSimulatorOutput(problemDG0);
[wsDG1, stDG1] = getPackedSimulatorOutput(problemDG1);
[wsWENO, stWENO] = getPackedSimulatorOutput(problemWENO);
[wsFIRef, stFIRef] = getPackedSimulatorOutput(problemFIRef);

%%

plotWellSols({wsDG0, wsDG1, wsWENO, wsFIRef})