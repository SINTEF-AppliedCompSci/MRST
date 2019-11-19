[wsSI, stSI]     = getPackedSimulatorOutput(problemSI);
[wsDG0, stDG0]     = getPackedSimulatorOutput(problemDG0);
[wsDG1, stDG1]     = getPackedSimulatorOutput(problemDG1);
[wsWENO, stWENO]   = getPackedSimulatorOutput(problemWENO);
[wsFIRef, stFIRef] = getPackedSimulatorOutput(problemFIRef);

%%

close all

G    = problemDG0.SimulatorSetup.model.G;
GRef = problemFIRef.SimulatorSetup.model.G;

step = 100;
figure
plotCellData(G, stDG0{step}.c)
axis equal tight
figure
plotCellData(G, stDG1{step}.c)
axis equal tight
figure
% plotCellData(G, stWENO{step}.c)
% figure
% plotCellData(GRef, stFIRef{step}.c)