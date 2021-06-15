mrstModule add ad-core ad-blackoil ad-props
%% Different options
% type = 'base';
% type = 'sparse'
type = 'diagonal';
useSparseBlocks = true;
useDiagonalMex = false;
useDiagonalOperators = true;
useDiagonalRowMajor = false;
useDiagonalDeferredAssembly = false;
if ~exist('n', 'var')
    n = 10;
end
p0 = 1*barsa;
s0 = [1, 1, 1]/3;
%%
G = cartGrid([n, n], [1, 1]);
G = computeGeometry(G);
rock = makeRock(G, 1*darcy, 0.5);
fluid = initSimpleADIFluid();
model = GenericBlackOilModel(G, rock, fluid);

switch type
    case 'sparse'
        backend = SparseAutoDiffBackend('useBlocks', useSparseBlocks);
    case 'diagonal'
        backend = DiagonalAutoDiffBackend('useMex',           useDiagonalMex, ...
                                      'rowMajor',         useDiagonalRowMajor, ...
                                      'deferredAssembly', useDiagonalDeferredAssembly, ...
                                      'modifyOperators',  useDiagonalOperators);
    otherwise
        backend = AutoDiffBackend();
end
model.AutoDiffBackend = backend;        % Set backend
model = model.validateModel();          % Validate model with new backend
state0 = initResSol(G, p0, s0);         % Set up initial state
state0 = model.validateState(state0);   % Validate initial state
state = model.getStateAD(state0);       % Initialize AD-state
forces = model.getValidDrivingForces(); % Set up dummy forces

%% Compute cached functions
timer = tic();
% Get two primary variables
[p, sw] = model.getProps(state, 'pressure', 'sw');
% Get the mass in each cell and the phase flux
[cm, v] = model.getProps(state, 'ComponentTotalMass', 'PhaseFlux');
eqs = model.getModelEquations(state0, state, 1*day, forces);
toc(timer)
%% Print statistics
fprintf('Grid has %d cells with %d interfaces\n', numelValue(p), numelValue(v{1}))
%% Show pressure variable with AD 
disp(p)
%% Show jacobian
Jp = p.jac{1}; Jw = sw.jac{1}; % Get pressure and S_w Jacobians
disp(Jp)                       % Show the pressure Jacobian
%% Pick some element
Jps = Jp(5, :) % Pick the fifth element of pressure Jacobian
%%
Jp(5, :) = 2*Jps; % Inserting in the same position is still diagonal
class(Jp)         % -> DiagonalJacobian: Diagonal structure preserved
Jp(6, :) = 3*Jps; % Insertion in another position results in sparse
class(Jp)         % -> double: Sparse matrix
%%
Jws = Jw(5, :);       % Take fifth element of S_w Jacobian
class(Jps + Jws)      % Different variables, same subset -> DiagonalJacobian
class(Jps + Jw(6, :)) % Different variables, different subset -> double
%%
Jp = sparse(Jp);
fprintf('Jacobian has type %s with dimensions %d by %d.\n', class(Jp), size(Jp))
%%
vw = v{1} % Get the water flux
%%
disp(vw.jac{1})
%%
ew = eqs{1} % Equation for conservation of water component
%%
disp(ew.jac{1})
%% Indicate why diagonal AD is a good idea
%% Memory usage
for i = 1:numel(p.jac)
    Jp = p.jac{i};
    whos J
end
%%
tic()
for i = 1:100
    tmp = p*5;
end
toc()
%%
if ~exist('m', 'var')
    m = 1e6;
end
if ~exist('nder', 'var')
    nder = 5;
end
its = 100;
Dx = rand(m, nder); Dy = rand(m, nder); % Diagonals of Jacobians
x = rand(m, 1); y = rand(m, 1);         % Vector values
fprintf('Testing with %d elements and %d derivatives per element\n', m, nder)

%% Sparse version
Jx = cell(1, nder); Jy = cell(1, nder);
for i = 1:nder
    Jx{i} = sparse(1:m, 1:m, Dx(:, i), m, m);  % m-by-m diagonal matrix, J_xi
    Jy{i} = sparse(1:m, 1:m, Dy(:, i), m, m);  % m-by-m diagonal matrix, J_yi
end
tic()
for it = 1:its
    J_sparse = cell(1, nder);
    dx = sparse(1:m, 1:m, x, m, m);
    dy = sparse(1:m, 1:m, x, m, m);
    for i = 1:nder
        J_sparse{i} = dx*Jy{i} + dy*Jx{i};
    end
end
t_sparse = toc()/it;
fprintf('Sparse: Average time is %1.4fs\n', t_sparse)

%% Diagonal version
tic()
for it = 1:its
    % Older Matlab versions
    % J_diag{i} = bsxfun(@times, D1, v2) + bsxfun(@times, D2, v1);
    J_diag = Dx.*y + Dy.*x;
end
t_diag = toc()/it;
fprintf('Diagonal: Average time is %1.4fs\n', t_diag)
%%
fprintf('Sparse: %1.4fs Diagonal: %1.4fs (%2.1f speed up from diagonal)\n', t_sparse, t_diag, t_sparse/t_diag)

%% Copyright Notice
%
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
