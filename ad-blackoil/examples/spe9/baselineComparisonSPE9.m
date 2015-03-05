%% Ninth Comparative Solution Project
% This example runs the model from Killough, J. E. 1995. Ninth SPE
% comparative solution project: A reexamination of black-oil simulation. In
% SPE Reservoir Simulation Symposium,  12-15 February 1995, San Antonio,
% Texas. SPE 29110-MS, doi: 10.2118/29110-MS

%% Set up model
% This <ninth SPE comparative solution project> consists of a water
% injection problem in a highly heterogenous reservoir. There is one
% injector and 25 producers. The problem is set up to be solved using a
% black-oil model. The data set we provide is a modified version of input
% files belonging to the <http://www.ntnu.edu/studies/courses/TPG4535
% course in reservoir engineering and petrophysics> at NTNU (Trondheim,
% Norway) and available at
% <http://www.ipt.ntnu.no/~kleppe/pub/SPE-COMPARATIVE/ECLIPSE_DATA/>.
%
% We have outsourced the setup of the schedule, model and initial state to
% the getBenchmarkAD function. 

mrstModule add ad-blackoil ad-core mrst-gui ad-unittest
[schedule, model, state0] = getBenchmarkAD('spe9');

%% Set up linear solver
% We proceed to setup a CPR-type
% solver, using the AGMG linear solver as the multigrid preconditioner. The
% CPR preconditioner attempts to decouple the fully implicit equation set
% into a pressure component and a transport component.
%
% The pressure is mathematically elliptic/parabolic in nature, and
% multigrid is well suited for solving these highly coupled, challenging
% problems. The remainder of the linear system, corresponding to the
% hyperbolic part of the equations is localized in nature and is primarily
% concerned with moving the saturation between neighboring grid blocks.
%
% Setting up a preconditioner is not strictly required to solve this
% problem, as the 9000 cells for a three-phase system results in linear
% systems with 27000 unknowns (one equation per phase, per cell). However,
% it does improve the solution speed and is required for larger cases,
% where Matlab's standard linear solvers scale poorly.
try
    mrstModule add agmg
    pressureSolver = AGMGSolverAD('tolerance', 1e-4);
catch
    pressureSolver = BackslashSolverAD();
end
linsolve = CPRSolverAD('ellipticSolver', pressureSolver);

%% Plot the rock permeability
% The SPE9 data set has a anisotropic, inhomogenous permeability field. The
% vertical permeability is 1/10th of the horizontal values. We plot the
% permeability using a log10 transform to better see the contrast.
v = [15, 30];
G = model.G;
rock = model.rock;

clf;
plotCellData(G, log10(rock.perm(:, 1)))
logColorbar();
axis tight
view(v)
title('SPE9 Horizontal permeability');

%% Plot the grid porosity
clf;
plotCellData(G, rock.poro)
colorbar();
axis tight
view(v)
title('SPE9 Porosity');
%% Plot a single vertical set of cells
% While the grid is structured, the grid has varying cell size along the
% vertical axis. To show this in detail, we plot the porosity in a single
% column of cells. We get the underlying logical grid and extract the
% subset corresponding to the first column (upper left corner of the grid).
%
% By using axis equal we see the actual aspect ratio of the cells.
clf;
[ii, jj] = gridLogicalIndices(G);
plotCellData(G, rock.poro, ii == 1 & jj == 1)
colorbar
axis equal tight
view(v);
%% Plot initial saturation of oil and water
% Initially, the reservoir does not contain free gas. We plot the initial
% saturations using the RGB plotting feature, where a three column matrix
% sent to plotCellData is interpreted column wise as fractions of red,
% green and blue respectively. Since MRST convention is to order the phases
% in the order WOG we need to permute the column index slightly to get red
% oil, blue water and green gas.
s = state0.s(:, [2, 3, 1]);
clf;
plotCellData(G, s)
axis tight
view(v)

%% Examine gas dissolved in oil phase
% Even though there is no free gas initially, there is significant amounts
% of gas dissolved in the oil phase. The dissolved gas will bubble into
% free gas when the pressure drops below the bubble point pressure. For a
% given pressure there is a fixed amount of gas which can be dissolved in
% the black-oil instantanous dissolution model. To illustrate how saturated
% the initial conditions are, we plot the function
%
% $$ g(p) = \frac{R_s}{R_s^{sat}(p)}  $$
%
% I.e. how close the oil phase is to being completely saturated for the
% current pressure. A value near 1 means that the liquid is close to
% saturated and any values above 1 will immediately lead to free gas
% appearing in the simulation model.
%
% As we can see from the figure, free gas will appear very quickly should
% the pressure drop.
Rs_sat = model.fluid.rsSat(state0.pressure);
Rs = state0.rs;

clf;
plotCellData(G, Rs./Rs_sat)
axis tight
colorbar
view(v)
title('Fraction of maximum gas saturation in oil phase - g(p)');

%% Run the schedule
[wellsols, states, reports] =...
    simulateScheduleAD(state0, model, schedule, 'LinearSolver', linsolve);

%% Load comparison data from commercial solver
addir = mrstPath('ad-blackoil');
compare = fullfile(addir, 'examples', 'spe9', 'compare');
[smry, smspec]  = readSummaryLocal(fullfile(compare, 'SPE9'));

compd = 1:(size(smry.data, 2));
Tcomp =  smry.get(':+:+:+:+', 'YEARS', compd);


%%
wsign = [wellsols{1}.sign];

inj = find(wsign == 1);
prod = find(wsign == -1);
T = convertTo(cumsum(schedule.step.val), year);

%%
mrstplot = @(data) plot(T, data, '-b', 'linewidth', 2);
compplot = @(data) plot(Tcomp, data, 'ro', 'linewidth', 2);
%%
clf
names = {'PROD13', 'PROD18'};
nn = numel(names);
for i = 1:nn

    name = names{i};
    
    comp = convertFrom(smry.get(name, 'WBHP', compd), psia)';
    mrst = getWellOutput(wellsols, 'bhp', name);
    
    subplot(nn, 1, i)
    hold on
    mrstplot(mrst);
    compplot(comp);
    title(name)
    axis tight
    grid on
    
    xlabel('Time (years)')
    ylabel('Pressure (Pa)')
end
legend({'MRST', 'ECL'})
%%
clf
for i = 1:nn
    name = names{i};
    comp = convertFrom(smry.get(name, 'WGPR', compd), 1000*ft^3/day);
    mrst = abs(getWellOutput(wellsols, 'qGs', name));
    
    subplot(nn, 1, i)
    hold on
    mrstplot(mrst);
    compplot(comp);
    title(name)
    axis tight
    grid on
    
    xlabel('Time (years)')
    ylabel('Gas rate (m^3/s)')
end
legend({'MRST', 'ECL'})

%% Changing controls

%% Saturation fronts etc
