%% Show Sparsity Pattern
% This example compares the sparsity patterns of dG(0)/SPU and dG(1) for
% a small PEBI grid and shows how one can use potential ordering to permute
% the systems to (block)triangular form.
mrstModule add upr dg ad-core ad-props ad-blackoil sequential

%% Set up the model and compute a pressure solution
% We construct a PEBI grid using a function from the UPR module, specify a
% quarter five-spot type well pattern, and set up an incompressible
% oil-water system with unit properties and slight rock compressibility.
% We then split the generic black-oil model into a pressure and a transport
% model (with dG) and use a standalone pressure solver to compute pressure.
% This is used as input to compute the linearized transport equations for
% the dG transport solver.
G        = pebiGrid2D(1/4, [1,1]);
G        = computeGeometry(G);
G        = computeCellDimensions(G);
rock     = makeRock(G, 1, 1);
W        = addWell([], G, rock, 1, 'type', 'rate', 'Radius', 1e-3, ...
                    'val',  1, 'compi', [1 0], 'Name', 'I');
W        = addWell(W, G, rock, G.cells.num, 'type', 'rate', 'Radius', 1e-3,...
                    'val', -1, 'compi', [1 0], 'Name', 'P');
fluid    = initSimpleADIFluid('phases', 'WO', 'mu', [1, 1], 'rho', [1, 1], 'cR', 1e-6/barsa);
schedule = simpleSchedule(1, 'W', W);

model    = GenericBlackOilModel(G, rock, fluid, 'gas', false);
pmodel   = PressureModel(model); % Pressure model
tmodel0  = TransportModelDG(model, 'degree', 0);
model0   = SequentialPressureTransportModel(pmodel, tmodel0,'parentModel', model);

forces   = schedule.control(1);
model0   = model0.validateModel(forces);
state0   = model0.validateState(initResSol(G, 1*barsa, [0, 1]));
state    = standaloneSolveAD(state0, model0, schedule.step.val(1), 'W', forces.W);
problem  = model0.transportModel.getEquations(state0, state, 1, forces);

%% Sparsity pattern for dG(0)/SPU
% Show sparsity patter before and after applying potential ordering
figure('Position',[520 360 985 420]);
subplot(1,2,1)
A0 = problem.equations{1}.jac{1};
spy(A0);
set(gca,'XTick',[],'YTick',[]); xlabel('');
axes('Position',[.31 .555 .15 .35]);
plotGrid(G,'FaceColor',[1 1 .9]);
h=text(G.cells.centroids(:,1),G.cells.centroids(:,2),num2str((1:G.cells.num)'));
set(h,'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off

subplot(1,2,2)
[~,q] = sort(state.pressure,'descend');
spy(A0(q,q));
set(gca,'XTick',[],'YTick',[]); xlabel('');
axes('Position',[.75 .555 .15 .35]);
plotGrid(G,'FaceColor',[1 1 .9]);
h=text(G.cells.centroids(q,1),G.cells.centroids(q,2),num2str((1:G.cells.num)'));
set(h,'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off

%% Get the discretization matrix for the dG(1) scheme
% As for the dG(0)/SPU scheme, we first construct a transport model and a
% corresponding sequential model, before we compute one pressure step and
% use the result as input to construct the linear transport system
tmodel1  = TransportModelDG(model, 'degree', 1);
model1   = SequentialPressureTransportModel(pmodel, tmodel1,'parentModel', model);

model1   = model1.validateModel(forces);
state0   = model1.validateState(initResSol(G, 1*barsa, [0, 1]));
state    = standaloneSolveAD(state0, model1, schedule.step.val(1), 'W', forces.W);
problem  = model1.transportModel.getEquations(state0, state, 1, forces);

%% Sparsity pattern for dG(1)
% Show the sparsity system before and after potential ordering
figure('Position',[520 360 985 420]);
subplot(1,2,1)
A1 = problem.equations{1}.jac{1};
spy(A1,8);
set(gca,'XTick',[],'YTick',[]); xlabel('');
axes('Position',[.31 .555 .15 .35]);
plotGrid(G,'FaceColor',[1 1 .9]);
h=text(G.cells.centroids(:,1),G.cells.centroids(:,2),num2str((1:G.cells.num)'));
set(h,'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off

subplot(1,2,2)
qq = 3*ones(G.cells.num,1); qq([1 end])=1; qq=cumsum([0; qq]);
map = cell(G.cells.num,1);
for i=1:G.cells.num 
    map{i} = (qq(i)+1:qq(i+1))';
end
p=vertcat(map{q});
spy(A1(p,p),8);
set(gca,'XTick',[],'YTick',[]); xlabel('');
axes('Position',[.75 .555 .15 .35]);
plotGrid(G,'FaceColor',[1 1 .9]);
h=text(G.cells.centroids(q,1),G.cells.centroids(q,2),num2str((1:G.cells.num)'));
set(h,'FontSize',8,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
