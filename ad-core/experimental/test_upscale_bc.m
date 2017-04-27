mrstModule add ad-props ad-core ad-blackoil
G = cartGrid([10, 10, 10]);
G = computeGeometry(G);

W = [];
src = [];
bc = [];
bc = pside(bc, G, 'xmax', 100*barsa, 'sat', [1, 0]);
bc = fluxside(bc, G, 'xmin', 100, 'sat', [1, 0]);

src = addSource(src, 1:G.cells.num, 1, 'sat', [1, 0]);
schedule = simpleSchedule(1, 'bc', bc, 'W', W, 'src', src);

p = partitionUI(G, [2, 2, 2]);
CG = generateCoarseGrid(G, p);

rock = makeRock(G, 1, 1);
fluid = initSimpleADIFluid();

model = TwoPhaseOilWaterModel(G, rock, fluid);
model_c = upscaleModelTPFA(model, p);

schedule_c = upscaleSchedule(model_c, schedule, 'bcUpscaleMethod', 'linear');

% src = [];
% addSource(