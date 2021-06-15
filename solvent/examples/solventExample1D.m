%% Using the black-oil solvent model
% In this example, we investigate black-oil solvent model, which models the
% effect of injecting a solvent gas that mixes with the reservoir
% hydrocarbon fluids. The model is an extension of the Todd-Longstaff
% mixing model, which treats the solvent gas as a fourth pseudo-component.
mrstModule add ad-core ad-props ad-blackoil spe10 mrst-gui solvent

df = get(0, 'defaultfigureposition');
close all

%% Set up grid and rock
% We will set up a simple 1D example with homogeneous rock properties
n = 50;
l = 100*meter;
G = computeGeometry(cartGrid([n,1,1], [l,1,1]));

perm = 100*milli*darcy;
poro = 0.4;
rock = makeRock(G, perm, poro);

%% Set up fluid and add solvent gas properties
% The solvent model we will use, treats solvent gas as a fourth
% pseudocomponent, which is either miscible or immiscible with the oil and
% gas, depending on the fraction of solvent concentration to total gas
% concentration, $S_s/(S_g + S_s)$, and the pressure. The model assumes one
% immiscible and one miscible residual saturation for the hydrocarbon
% phases, with $S_{\alpha r,i} > S_{\alpha r,m}$, and we must define both.
% In this example, we assume that the critical (residual) gas saturation is
% zero in both the immiscible and miscible case. The degree of mixing is
% modeled through a mixing paramter $\omega$ that defines the degree of
% mixing, where (no mixing) $= 0 \leq \omega \leq 1 =$ (full mixing).

% Set up three-phase fluid with quadratic relperms
fluid = initSimpleADIFluid('n'     , [2, 2, 2]                        , ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG'                            , ...
                           'mu'    , [1, 3, 0.4]*centi*poise          );

sOr_i = 0.30; % Immiscible residual oil saturation
sOr_m = 0.10; % Miscible residual oil saturation
sGc_i = 0.15; % Immiscible residual oil saturation
sGc_m = 0.05; % Miscible residual oil saturation

% We scale the relperms to the immiscible endpoints
fluid.krW = coreyPhaseRelpermAD(2, 0    , fluid.krW(1 - (sOr_i + sGc_i)), sOr_i + sGc_i);
fluid.krG = coreyPhaseRelpermAD(2, sGc_i, fluid.krG(1 - sOr_i), sOr_i + sGc_i);     
[fluid.krO, fluid.krOW, fluid.krOG] ...
     = deal(coreyPhaseRelpermAD(2, sOr_i, fluid.krO(1 - sGc_i), sOr_i + sGc_i));

% Add the solvent pseudocomponent to the fluid
fluid   = addSolventProperties(fluid, 'rhoSS' , 90*kilogram/meter^3, ...
                                      'muS'   , 0.5*centi*poise    , ...
                                      'mixPar', 1                  , ...
                                      'sOr_m' , sOr_m              , ...
                                      'sGc_m' , sGc_m              );
                                  
%% Set up black-oil solvent model
% Since the residual saturations changes with the solvent and pressure,
% dynamic endpoint-scaling should be used in the relperm evaluations.
% However, this is somewhate unstable in cases where the initial oil or gas
% saturation is close to the residual saturation, and is thus disabeled by
% default.
model_fm = BlackOilSolventModel(G, rock, fluid                   , ...
                                       'extraStateOutput'      , true, ...
                                       'dynamicEndPointScaling', true);

%% Inspect fluid model
% In regions with only oil, reservoir gas and water, we have traditional
% black-oil behavior, whereas regions with only water, reservoir oil and
% solvent gas, the oil is completely miscible with the solvent. In this
% case, solvent gas mixes with formation oil, effectively altering the
% viscosities, densities and relperms according the the Todd-Longstaff
% model [1]. In the intermediate region, we interpolate between the two
% extrema depending on the solvent to total gas saturation fraction, and
% the pressure. We look at terneary plots for the hydrocarbon relperms and
% viscosities when no water is present. Notice how the residual oil
% saturation (white line) reduces from $S_{or,i}$ (immiscible) to
% $S_{or,m}$ (miscible) with increasing solvent saturation, and the kinks
% in the gas and solvent relperms across this line.
plotSolventFluidProps(model_fm, {'kr', 'mu'}, {'O', 'G', 'S'});

%% Set up wells and simulate schedule
% We set up a simple injection schedule where solvent is injected into the
% reservoir at constant rate over a period of three years. The reservoir
% is initially filled with a mixture of oil and water.
time = 3*year;
dt   = 30*day;
rate = 1*sum(poreVolume(G, rock))/time;
bhp  = 50*barsa;
W    = [];
W    = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate, 'comp_i', [0,0,0,1]);
W    = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', bhp , 'comp_i', [1,0,0,0]);

tvec     = rampupTimesteps(time, dt);
schedule = simpleSchedule(tvec, 'W', W);

state0         = initResSol(G, 100*barsa, [0.6, 0.4, 0, 0]);
state0.wellSol = initWellSolAD(W, model_fm, state0);
[wellSols_fm, states_fm, reports_fm] = simulateScheduleAD(state0, model_fm, schedule);

%% Compare with moderate mixing
% To see the effect of the mixing parameter $\omega$, we compare the
% results above with a model with moderate mixing, where we set $\omega =
% 1/3$.
model_mm = model_fm;
model_mm.fluid.mixPar = 1/3;
[wellSols_mm, states_mm, reports_mm] = simulateScheduleAD(state0, model_mm, schedule);

%% Plot the results
% We now compare the results obtained from the two models by plotting the
% oil and solvent mass in each cell for three time steps. We see how the
% model with full mixing predicts a much higher oil recovery, and that the
% solvent travels faster with higher mixing due to the lower viscosity and
% thus higher mobility of the solvent-oil mixture. The difference in oil
% mass is shown as grey-shaded areas.

% Plotting parameters
x      = linspace(0,l,n)';
phases = [2,4];
clr    = lines(numel(phases));
mrk    = {'-', '--'};
pargs  = @(phNo, mNo) {mrk{mNo}, 'color', clr(phNo,:), 'lineWidth', 2};
step   = [15, 25, 45];
pv     = model_fm.operators.pv;
gray   = [1,1,1]*0.9;
time   = cumsum(schedule.step.val);

% Residual oil mass in the immiscible and miscible case
massr_i = pv.*sOr_i*fluid.rhoOS;
massr_m = pv.*sOr_m*fluid.rhoSS;

% Mass of solvent and gas phase for each time step
getMass = @(states) cellfun(@(state) bsxfun(@times, pv, state.s(:,phases).*state.rho(:,phases)), ...
                                           states, 'uniformOutput', false);
mass = {getMass(states_fm); getMass(states_mm)};

figure('Position', [df(1:2), 1000, 400]);

for sNo = 1:numel(step)
    subplot(1, numel(step),sNo);
    box on; hold on
    for phNo = 1:numel(phases)
        if phases(phNo) == 2
                xx = [x; x(end:-1:1)];
                mm = [mass{1}{step(sNo)}(:, phNo); 
                      mass{2}{step(sNo)}(end:-1:1, phNo)];
                fill(xx, mm, gray, 'edgecolor', 'none');
        end
        for mNo = 2:-1:1
            prg = pargs(phNo, mNo);
            [hm(mNo), hp(phNo)] = deal(plot(x, mass{mNo}{step(sNo)}(:,phNo), prg{:}));
        end
    end
    hold off
    if sNo == 1
        legend(hp, {'Oil mass', 'Solvent mass'}, 'Location', 'NorthWest')
        ylabel('Mass [kg]');
    elseif sNo == 2
        xlabel('Distance from injector [m]');
    elseif sNo == numel(step)
        legend(hm, {'Full mixing', 'Moderate mixing'}, 'Location', 'NorthWest')
    end
    ylim([0,650]);
    title([num2str(time(step(sNo))/day), ' days'])
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}