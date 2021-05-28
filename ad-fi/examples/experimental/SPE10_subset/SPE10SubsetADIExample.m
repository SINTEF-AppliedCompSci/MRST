%% Simulate a large example using parts of SPE10
% This example is a larger example demonstrating the solver on a medium
% size grid (66,000 cells) with a relatively large amount of time steps
% (100). This example will take some time, especially if MLDIVIDE is used
% as the elliptic solver. Be wary that increasing the number of layers may
% let the simulations take a very long time.
%
% Furthermore, note that the 'solvefiADI' solver and especially the CPR
% preconditioner, uses direct indexing into sparse matrices.  In MATLABs
% prior to release R2011a, this is relatively inefficient.  Starting from
% 2011a however, the bottleneck of direct indexing has been largely
% removed.
mrstModule add ad-core ad-fi deckformat spe10

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE10-S3.DATA.txt');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);


% Create a special ADI fluid which can produce differentiated fluid
% properties.
fluid = initDeckADIFluid(deck);

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
fluid.krO = fluid.krOW;

% The case includes gravity
gravity on


% The initial state is provided as a binary file. The initial state
% contains a uniform mixture of water (.12) and oil (.88).
% load initialState;
%% Set up permeability, grid and wells
% We will simulate on the top 5 layers.
layers = 1:5;

[G, W, rock] = SPE10_setup(layers);


% SPE10 contains zero and extremely low porosities. For the purpose of this
% tutorial, we will mask away these values. An alternative would be to set
% these cells to inactive by using extractSubgrid and removing the
% corresponding cells.
low = 1e-4;
rock.poro(rock.poro < low) = low;

%% Plot the permeability
clf;
plotCellData(G, log10(rock.perm(:,1)));

%% Set up solution structures.

% The initial reservoir is at 6000 psi and is fully oil saturated. The well
% solution gets its initial pressure from the bottom hole pressure values
% provided.
initSat = [0 1];
state0 = initResSol(G, 6000*psia, initSat);
state0.wellSol = initWellSolLocal(W, state0);

for i = 1:numel(W)
    state0.wellSol(i).bhp = W(i).val;
    % Set well sign
    if strcmpi(W(i).name(1), 'p')
        W(i).sign = -1;
    else
        W(i).sign = 1;
    end
end

% Set up a Water / Oil system using CPR preconditioner. Alternatively we
% could have used a specialized elliptic solver using the option 'cprEllipticSolver'
% to exploit the nature of the variables involving pressure.

system = initADISystem({'Water', 'Oil'}, G, rock, fluid, 'cpr', true, 'cprRelTol', 2e-2);
% If an alternative solver for the pressure subproblem was installed, it
% could be added using
%
%    system.nonlinear.cprEllipticSolver = @(A,b) solver(A,b)
%
% This can greatly speed up the solution process, as the default option
% uses MATLABs direct solver @mldivide which is very expensive for a
% preconditioner.

%% Simulate 1000 days of production and save iteration count and time
% We provide the solver with time steps for roughly 1000 days of
% production. A few smaller steps are done to get better accuracy during
% the initial injection period. After this we do 10 day intervals to step
% rapidly through the schedule. While this converges at every time step,
% implicit solvers will still get improved accuracy by doing smaller time
% steps. Numerical diffusion can, for instance, be problematic when doing
% large time steps.

dt = [        0.001*day;
      repmat( 0.1  *day, [  5, 1]);
      repmat( 1    *day, [ 10, 1]);
      repmat(10    *day, [100, 1])];
nstep = numel(dt);

states = cell(nstep,1);
its = zeros(nstep,1);
time = zeros(nstep,1);

state = state0;
for t = 1 : nstep
    fprintf('Step %d/%d: %.2f -> %.2f [days]\n\n', t, nstep, ...
            convertTo(sum(dt(1:t-1)), day), ...
            convertTo(sum(dt(1:t  )), day));

    timer = tic();

    [state, it] = solvefiADI(state, dt(t), W, G, system);

    states{t} = state;
    its(t) = it;
    time(t) = toc(timer);

    fprintf('\n\n');
end

%%
h = [];
set(gcf, 'Renderer', 'OpenGL');

[htop, htext, hs] = plotWell(G, W);
hg = plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .05);
view(-60,  70);
axis tight off

%% Plot the solution
for i = 1:numel(states)
    delete(h);
    data = states{i}.s(:,1);
    h = plotGrid(G, data > 0, 'FaceAlpha', .3, ...
                 'FaceColor', 'red', 'EdgeAlpha', 0);
    % Uncomment for slower, but prettier volume plotting
%     plotGridVolumes(G, data);
    title(['Water front after ' formatTimeRange(sum(dt(1:i)))])
    pause(.1)
end

%% Plot the time taken and number of iterations

delete([h, htext, htop, hs, hg])

plot([time, its], '.-');
legend({'Time (s)', 'Iterations'})
xlabel('Step #')
title('Time and iteration count')

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
