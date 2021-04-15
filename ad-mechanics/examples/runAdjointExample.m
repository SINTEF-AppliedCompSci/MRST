%% Computation of Gradients using Adjoint simulations
%
% In this example, we demonstrate how one can setup adjoint simulations in a
% poroelastic simulation in order to compute gradients (or sensitivities) of a
% given quantity. Here, we will compute the gradient of the vertical
% displacement (uplift) at a node with respect to the injection rates.
%
% We consider a 2D domain, with two injection wells located on the sides and one
% production well in the middle.
%
% The injection rate is given by a schedule (see below). At the production
% well, we impose a constant pressure.
%
% We compute the gradient of the uplift (vertical displacement) at the top of
% the domain with respect to the injection rates at every time step.
%
% We use this gradient information to update the schedule so that the uplift is
% reduced. Different measures of the time average for the uplift are considered,
% and they produce different result, see the last plot in this script.
%
% Different fluid model can be used (two phases or single phase).
%
%

mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Setup geometry
%

% We consider a  2D regular cartesian domain

cartDim = [31, 30];
L       = [30, 10];
G = cartGrid(cartDim, L);
G = computeGeometry(G);

%% Setup fluid

% We can consider different fluid models

opt.fluid_model = 'single phase';
pRef = 100*barsa;
switch opt.fluid_model
  case 'single phase'
    fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                               1000*kilogram/meter^3, 'c', 1e-4, 'cR', ...
                               4e-10, 'pRef', pRef);
  case 'oil water'
    fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1, 100]*centi*poise, 'n', ...
                               [1, 1], 'rho', [1000, 700]*kilogram/ meter^2, 'c', ...
                               1e-10*[1, 1], 'cR', 1e-10, 'pRef', pRef);
  case 'blackoil'
    error('not yet implemented, but could be easily done!')
  otherwise
    error('fluid_model not recognized.')
end



%% Setup rock parameters (for flow)
%

rock.perm = darcy*ones(G.cells.num, 1);
rock.poro = 0.3*ones(G.cells.num, 1);


%% Setup material parameters for Biot and mechanics
%

E          = 1e-2*giga*Pascal; % Young's module
nu         = 0.3;              % Poisson's ratio
alpha      = 1;                % Biot's coefficient
% Convert global properties to cell values
E          = repmat(E, G.cells.num, 1);
nu         = repmat(nu, G.cells.num, 1);
rock.alpha = repmat(alpha, G.cells.num, 1);


%% Setup boundary conditions for mechanics (no displacement)
%
%
% zero displacement at bottom, left and right sides. We impose a given pressure
% at the top.

% Gather the Dirichlet boundary faces (zero displacement) at left, bottom and right.
dummyval = 100; % We use pside to recover face at bottom, we use a dummy
                % value for pressure in this function.
bc = pside([], G, 'Xmin', dummyval);
bc = pside(bc, G, 'Xmax', dummyval);
bc = pside(bc, G, 'Ymin', dummyval);
indfacebc = bc.face;

% Get the nodes that belong to the Dirichlet boundary faces.
facetonode = accumarray([G.faces.nodes, rldecode((1 : G.faces.num)', ...
                                                 diff(G.faces.nodePos))], ...
                        ones(numel(G.faces.nodes), 1), [G.nodes.num, ...
                    G.faces.num]);
isbcface = zeros(G.faces.num, 1);
isbcface(indfacebc) = 1;
bcnodes  = find(facetonode*isbcface);
nn       = numel(bcnodes);
u        = zeros(nn, G.griddim);
m        = ones(nn,  G.griddim);
disp_bc  = struct('nodes', bcnodes, 'uu', u, 'mask', m);

% Set a given pressure on the  face at the top.
dummyval = 100; % We use pside to recover face at bottom, we use a dummy
                % value for pressure in this function.
bc = pside([], G, 'Ymax', dummyval);
sidefaces = bc.face;
signcoef = (G.faces.neighbors(sidefaces, 1) == 0) - (G.faces.neighbors(sidefaces, ...
                                                  2) == 0);
n = bsxfun(@times, G.faces.normals(sidefaces, :), signcoef./ ...
           G.faces.areas(sidefaces));
force = bsxfun(@times, n, pRef);
force_bc = struct('faces', sidefaces, 'force', force);

% Construct the boundary conidtion structure for the  mechanical system
el_bc = struct('disp_bc' , disp_bc, 'force_bc', force_bc);


%% Setup volumetric load for mechanics
%
% In this example we do not impose any volumetric force
loadfun = @(x) (0*x);


%% Gather all the mechanical parameters in a struct
%

mech = struct('E', E, 'nu', nu, 'el_bc', el_bc, 'load', loadfun);


%% Set gravity off
%

gravity off


%% Setup model
%

switch opt.fluid_model
  case 'single phase'
    model = MechWaterModel(G, rock, fluid, mech);
  case 'oil water'
    model = MechOilWaterModel(G, rock, fluid, mech);
  case 'blackoil'
    error('not yet implemented')
  otherwise
    error('fluid_model not recognized.')
end


%% Set up initial reservoir state
%
% The initial fluid pressure is set to a constant. We have zero displacement
% initially.
%

clear initState;
initState.pressure = pRef*ones(G.cells.num, 1);
switch opt.fluid_model
  case 'single phase'
    init_sat = [1];
  case 'oil water'
    init_sat = [0, 1];
  case 'blackoil'
    error('not yet implemented')
    % init_sat = [0, 1, 0];
    % initState.rs  = 0.5*fluid.rsSat(initState.pressure);
  otherwise
    error('fluid_model not recognized.')
end
% set up initial saturations
initState.s  = ones(G.cells.num, 1)*init_sat;
initState.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
% We compute the corresponding displacement field using the dedicated
% function computeInitDisp (actually not necessary as the solution should be zero).
initState    = computeInitDisp(model, initState, [], 'pressure', initState.pressure);
initState    = addDerivedQuantities(model.mechModel, initState);


%% Setup the wells
%
% Two injection wells on the sides and one production well in the middle.
%

nx = G.cartDims(1);
ny = G.cartDims(2);
switch opt.fluid_model
  case 'single phase'
    comp_inj  = [1];
    comp_prod = [1];
  case 'oil water'
    comp_inj  = [1, 0];
    comp_prod = [0, 1];
  case 'blackoil'
    error('not yet implemented')
  otherwise
    error('fluid_model not recognized.')
end

W = [];
wellopt = {'type', 'rate', 'Sign', 1, 'comp_i', comp_inj};
% Two injection wells vertically aligned, near the bottom
W = addWell(W, G, rock, round(nx/4)   + floor(1/4*ny)*nx, wellopt{:});
W = addWell(W, G, rock, nx + 1 - round(nx/4) + floor(1/4*ny)*nx, wellopt{:});
% production well in the center
wellopt = {'type', 'bhp', 'val', pRef, 'Sign', -1, 'comp_i', comp_prod};
W = addWell(W, G, rock, round(nx/2)   + floor(1/4*ny)*nx, wellopt{:});

% We plot the well location
wellcells = zeros(G.cells.num, 1);
wellcells(W(1).cells) = 1;
wellcells(W(2).cells) = 1;
wellcells(W(3).cells) = 2;
figure
clf
plotCellData(G, wellcells);
comment = ['The connection of the wells are colored' char(10) 'We have two ' ...
          'injection wells and on production well'];
text(0.5, 0.9, comment, 'units', 'normalized', 'horizontalalignment', 'center', ...
     'backgroundcolor', 'white');

% We incorporate the well in a FacilityModel which takes care of coupling all
% the well equations with the reservoir equations

facilityModel = FacilityModel(model.fluidModel);
facilityModel = facilityModel.setupWells(W);
model.FacilityModel = facilityModel;
model = model.validateModel(); % setup consistent fields for model (in
                               % particular the facility model for the fluid
                               % submodel)

%% Setup a schedule
%
%
% We set up a schedule where we gradually decrease from a maximum to a minimum
% injection rate value. Then, we keep the injection rate constant. See plot below.
%

clear schedule
schedule.step.val = [1*day*ones(1, 1); 10*day*ones(30, 1)];
nsteps = numel(schedule.step.val);
schedule.step.control = (1 : nsteps)';
valmax = 1*meter^3/day;
valmin = 1e-1*meter^3/day;
ctime = cumsum(schedule.step.val);
flattentime = 150*day; % Time when we reach the minimal rate value
qW = zeros(nsteps, 1);
for i = 1 : numel(schedule.step.control)
    if ctime(i) < flattentime
        qW(i) = valmin*ctime(i)/flattentime + valmax*(flattentime - ctime(i))/flattentime;
    else
        qW(i) = valmin;
    end
    W(1).val = qW(i);
    W(2).val = qW(i);
    schedule.control(i) = struct('W', W);
end
% We plot the injection schedule
figure
plot(ctime/day, qW*day);
axis([0, ctime(end)/day, 0, 1])
title('Initial schedule');
xlabel('time (day)')
ylabel('Injection rate (m^3/day)');


%% Run the schedule
%
% We run the simulation for the given model, initial state and schedule.
%

[wellSols, states] = simulateScheduleAD(initState, model, schedule);

% We start a visualization tool to inspect the result of the simulation
figure
plotToolbar(G, states);
colorbar


%% We plot the evolution of the uplift
%

% Get index of a node belonging to the cell at the middle on the top layer.
topcell = floor(nx/2) + nx*(ny - 1);
topface = G.cells.faces(G.cells.facePos(topcell) : (G.cells.facePos(topcell + ...
                                                  1)  - 1), :);
topface = topface(topface(:, 2) == 4, 1);
topnode = G.faces.nodes(G.faces.nodePos(topface)); % takes one node from the
                                                   % top face (the first listed)

laststep = numel(states);
uplift = @(step) (computeUpliftForState(model, states{step}, topnode));
uplifts = arrayfun(uplift, (1 : laststep)');
figure
plot(ctime/day, uplifts, 'o-');
title('uplift values');
xlabel('time (in days)');


%% Setup of the objective function given by averaged uplift values
%
% The adoint framework is general and computes the derivative of a (scalar)
% objective function of following forms,
%
%       obj = sum_{time step i = 1,..., last} (obj_i(state at time step i))
%
% This special form of this objective function allows for a recursive computation
% of the adjoint variables.
%
% We set up our objective function as being a time average of weighted uplift in
% a node in the middle at the top
%
% Given an exponent p, we take
%
%      obj = sum_{time step i = 1, ..., las} ( dt(i) * uplift(i) ^ p )
%
% See function objUpliftAver.
%
% We consider three values for the exponent p.


%% Compute gradients using the adjoint formulation for three different exponents


exponents = {1; 10};
np = numel(exponents);
adjointGradients = cell(np, 1);

for p = 1 : np
    C = 1e3; % For large exponents we need to scale the values to avoid very large
             % or very small objective function
    objUpliftFunc = @(tstep, ~, state) objUpliftAver(model, [], schedule, topnode, 'computeAllSteps', false, 'tStep', ...
                                                     tstep, 'state', state, 'computePartials', true, 'exponent', ...
                                                     exponents{p}, 'normalizationConstant', C);
    fprintf('\n***\n*  Start adjoint simulation for exponent p=%g\n*\n', exponents{p});
    adjointGradients{p} = computeGradientAdjointAD(initState, states, model, schedule, objUpliftFunc);
end

%% We can check the results from the adjoint computation by using finite difference
%
%
% The function computeGradientAdjointAD sets up this computation check for us.
% It should be used with a smaller schedule, otherwise the computation is very
% long.

compute_numerical_derivative = false;
if compute_numerical_derivative
    p = 3;
    exponent = exponents{p};
    objUpliftFunc2 = @(wellSols, states, schedule) objUpliftAver(model, states, schedule, topnode, 'computePartials', ...
                                                      false, 'exponent', exponent, 'normalizationConstant', C);
    fdgrad = computeGradientPerturbationAD(initState, model, schedule, objUpliftFunc2, 'perturbation', 1e-7);
    
    adjgrad = cell2mat(adjointGradients{p});
    fdgrad = cell2mat(fdgrad);
    
    abs(adjgrad - fdgrad)./(abs(adjgrad) + abs(fdgrad))

end

%% Plots of the results.
%
% We plot the uplift values with the gradients (normalized for comparison).

figure
clf
plot(ctime/day, uplifts, 'o-');
ylabel('uplift values');
xlabel('time (in days)');
hold on
yyaxis right

legendtext = {'uplift value'};
set(gca,'ColorOrder',hsv(3));
for p = 1 : np
    grads = cell2mat(adjointGradients{p});
    qWgrad = grads(1, :);
    % we renormalize the gradients to compare the two series of values
    qWgrad = 1/max(qWgrad)*qWgrad;
    plot(ctime/day, qWgrad, '*-');
    legendtext{end + 1} = sprintf('p = %g', exponents{p});
end

ylabel('gradient value')
legend(legendtext);


%% Update of the schedule based on gradient values
%
% We use the gradient values computed above to update the schedule so that the
% uplift, as measured with the three different values of the exponents, is
% decreased. The absolute minimum uplift is of course zero and it is obtained
% when we do not inject anything. To avoid this trivial case, we impose the
% extra requirement that, in all cases, we will have the same decrease in the total
% injection rates (see dcqW).

inituplifts = uplifts;
initschedule = schedule;
dts = schedule.step.val;
nsteps = numel(dts);
ctime = cumsum(schedule.step.val);
tottime = ctime(end);

np      = numel(adjointGradients);
uplifts = cell(np, 1);
qWs     = cell(np, 1);

dcqW = 3e-2*meter^3/day*tottime; % this is the total amount of injected fluid we
                                 % are willing to decrease.

for p = 1 : np

    grads = cell2mat(adjointGradients{p});

    qW = zeros(nsteps, 1);
    dqW = cell(2, 1);
    for i = 1 : 2
        dqW{i} = grads(i, :);
        dqW{i} = dcqW/(dqW{i}*dts)*dqW{i};
    end

    schedule = initschedule;
    for i = 1 : numel(schedule.step.control)
        W = schedule.control(i).W;
        % we update the values of the injection rate at each time step using
        % the gradient values.
        W(1).val = W(1).val - dqW{1}(i);
        W(2).val = W(2).val - dqW{2}(i);
        schedule.control(i) = struct('W', W);
        qW(i) = W(1).val;
    end

    qWs{p} = qW;

    figure
    clf
    plot(ctime/day, qW*day);
    title(sprintf(' Schedule obtained using p=%g', exponents{p}));
    xlabel('time (day)');
    ylabel('Injection rate (m^3/day)');
    
    fprintf('\n***\n*  Start simulation and compute uplift for updated schedule using exponent p=%g\n*\n', exponents{p});
    [wellSols, states] = simulateScheduleAD(initState, model, schedule);

    nsteps = numel(states);
    uplift = @(step) (computeUpliftForState(model, states{step}, topnode));
    uplifts{p} = arrayfun(uplift, (1 : nsteps)');

end

%% plot of the results.
figure
set(gcf, 'position', [100, 100, 1500, 600]);
clf
subplot(1, 2, 2)
hold on
plot(ctime/day, inituplifts);
legendtext = {'original uplift value'};
for p = 1 : np
    plot(ctime/day, uplifts{p});
    legendtext{end + 1} = sprintf('updated uplift (p=%g)', exponents{p});
end
xlabel('time (in days)');
legend(legendtext);

% Add some comments to the plot
[y, ind] = max(inituplifts);
x        = ctime(ind)/day;
set(gca, 'position', [0.5, 0.1, 0.4, 0.8], 'units', 'normalized')
xlim     = get(gca, 'xlim');
ylim     = get(gca, 'ylim');
apos  = get(gca, 'position');
xn = apos(1) + (x - xlim(1))/(xlim(2) - xlim(1))*apos(3);
yn = apos(2) + (y - ylim(1))/(ylim(2) - ylim(1))*apos(4);
annotation('textarrow', [xn, xn], [0.5, yn], 'string', 'maximum uplift');

subplot(1, 2, 1)
axis off
comments =  fileread('commentsToPlot.tex');
annotation('textbox', 'units', 'normalized', 'position', [0.05 0.1 0.4 0.8], ...
            'interpreter', 'tex', 'string', comments)

%%
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
