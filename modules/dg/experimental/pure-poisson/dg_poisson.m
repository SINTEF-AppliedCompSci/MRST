clear all
close all

mrstModule add vem dg vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection

mrstVerbose on;

% Basis polynomial degree
degree = 1;

% Symmetric or non-symmetric discretization
alpha = -1; % symmetric
%alpha = 1; % nonsymmetric (nitsche_param can be 0)

nitsche_param = 10;
if alpha == 1
  nitsche_param = 0;
end

% Number of elements
N = 10;

% Domain dimensions
L = 100;
h = L/N;

gridtype = "quad";
%gridtype = "triangle";
%gridtype = "pebi";

switch gridtype
  case "pebi"
    G = pebiGrid(h, [L, L]);
  case "triangle"
    G = cartGrid([N, N], [L, L]);
    G = triangleGrid(G.nodes.coords);
  case "quad"
    G = cartGrid([N, N], [L, L]);
  otherwise
    error(["unknown gridtype ", gridtype])
end
G = computeVEMGeometry(G);
G = computeCellDimensions(G);
%plotGrid(G),pause


% Rock porosity and permeability
rock = makeRock(G, 100*milli*darcy, 1);

% Fluid params
fluid = initSimpleADIFluid('phases', 'WO', ...
			   'rho'   , [1000, 1]*kilogram/meter^3, ...
			   'mu'    , [0.5, 0.5]*centi*poise , ...
			   'n'     , [1, 1]);


% Fluid model
modelfi = TwoPhaseOilWaterModel(G, rock, fluid);

% Transport and pressure models
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

% DG mrst discretization stuff
[jt, ot, mt] = deal(Inf);
jt = 0.2;
ot = 0.1;
mt = 0;
disc = DGDiscretization(modelDG.transportModel, ...
        'degree', degree, ...
        'basis', 'legendre', ...
        'useUnstructCubature', false, ... % Compute QR points with appropriate weights
        'jumpTolerance', jt, 'outTolerance', ot, 'meanTolerance', mt);

% Utilities
issymm = @(A) abs(max(max(A-A'))) < 1e-10;
equal = @(A, B) abs(max(max(A-B))) < 1e-10;

% Dofs and basis convenience functions
ndof = disc.basis.nDof;
cell_dofs = @(e) mcolon((e-1)*ndof+1, e*ndof)';
psi = disc.basis.psi;
gradpsi = disc.basis.grad_psi;

% Allocate scalar problem
ndofs = G.cells.num*ndof

% (grad v, grad w)
A = sparse(ndofs, ndofs);

% (<n.grad v>, [w])
S = sparse(ndofs, ndofs);

% Nitsche penalty
nitsche_penalty = @(h) nitsche_param*degree*degree/h;
P = sparse(ndofs, ndofs);

% Rhs
b = zeros(ndofs, 1);
loadfcn = @(x) 2*(pi/L)^2*sin(pi*x(:,1)/L).*sin(pi*x(:,2)/L);

% Fix exterior normals
normal_orientation_sign = @(G, e, f) sign(dot(mean(G.nodes.coords(G.faces.nodes(G.faces.nodePos(f):G.faces.nodePos(f+1)-1, :), :)) - mean(G.nodes.coords(G.cells.nodes(G.cells.nodePos(e):G.cells.nodePos(e+1)-1,:), :)), G.faces.normals(f,:)));
normal = @(G, e, f) normal_orientation_sign(G, e, f)*G.faces.normals(f,:) / norm(G.faces.normals(f,:));

% Test quadrature
exterior_area = 0;
interior_area = 0;

% Assemble
disp('Assemble...')
for e = 1:G.cells.num
  %disp(e)

  % Dofs
  dofs = cell_dofs(e);

  % Volume qr
  [Wvol, xvol] = disc.getCubature(e, 'volume');
  Wvol = full(Wvol)';
  nq = size(Wvol, 1);
  assert(abs(sum(Wvol)-G.cells.volumes(e)) < 1e-10)

  % Map to [-1,1]^2
  [xhatvol,~,scalevol] = disc.transformCoords(xvol, e);
  scalevolR = repmat(scalevol, nq, 1);

  % Get gradients
  gradpsi_i = zeros(nq, 2*ndof);
  for i = 1:ndof
    gradpsi_i(:, [2*i-1, 2*i]) = gradpsi{i}(xhatvol).*scalevolR;
  end

  % Put derivatives in same column
  gradpsi_i = reshape(gradpsi_i, 2*nq, ndof);

  % Fix size for .*
  WvolR = repmat(Wvol, 2, ndof);

  % Form stiffness
  A(dofs, dofs) = A(dofs, dofs) + gradpsi_i'*(gradpsi_i.*WvolR);

  % Local load
  psi_w = zeros(nq, ndof);
  for i = 1:ndof
    psi_w(:,i) = psi{i}(xhatvol).*Wvol;
  end
  b(dofs) = b(dofs) + psi_w'*loadfcn(xvol);

  % Loop over faces
  faces = G.cells.faces(G.cells.facePos(e):G.cells.facePos(e+1)-1, :);
  for f = faces(:,1)'
    cells = G.faces.neighbors(f, :);
    cell_idx = cells > 0;
    ncells = sum(cell_idx);
    cells = cells(cell_idx);

    % If no neighbor, let cells = [e e]
    if ncells == 1
      cells(2) = e;
    else
      % Make cells(1) == e
      if cells(1) ~= e
	cells = flip(cells);
      end
    end
    dofs = cell_dofs(cells);

    % Face QR
    [Wface, xface] = disc.getCubature(f, 'face');
    Wface = full(Wface)';
    nq = size(Wface, 1);
    assert(abs(sum(Wface) - G.faces.areas(f)) < 1e-10);
    [xhatface1,~,scale1] = disc.transformCoords(xface, cells(1));
    [xhatface2,~,scale2] = disc.transformCoords(xface, cells(2));

    % Debugging
    if ncells == 1
      exterior_area = exterior_area + sum(Wface);
    else
      interior_area = interior_area + 0.5*sum(Wface);
    end

    % Scale
    scale1R = repmat(scale1, nq, 1);
    scale2R = repmat(scale2, nq, 1);

    % Jumps and averages
    % NB: Only /2 if it's a shared f. Otherwise it's a one sided average.
    n = repmat(normal(G, e, f), nq, 1);
    jumps = zeros(nq, 2*ndof);
    averages = zeros(nq, 2*ndof);
    for i = 1:ndof
      jumps(:, [i ndof+i]) = [psi{i}(xhatface1) -psi{i}(xhatface2)];
      averages(:,[i ndof+i]) = [dot(gradpsi{i}(xhatface1).*scale1R, n, 2) ...
				   dot(gradpsi{i}(xhatface2).*scale2R, n, 2)]/ncells;
    end

    % Fix size for .*
    WfaceR = repmat(Wface, 1, 2*ndof);

    % -(<n.grad u>,[v])_F (possibly *2 if not a shared face (no averaging))
    SE = averages'*(jumps.*WfaceR);
    %SE = factor*jumps'*(averages.*WfaceR);

    % (beta/h [u],[v])_F
    PE = nitsche_penalty(G.cells.diameters(e))*jumps.'*(jumps.*WfaceR);

    % Assemble (possibly /2 for a shared face)
    locdofs = 1:ncells*ndof;
    dofs = dofs(locdofs);
    S(dofs, dofs) = S(dofs, dofs) + SE(locdofs, locdofs)/ncells;
    P(dofs, dofs) = P(dofs, dofs) + PE(locdofs, locdofs)/ncells;
  end

end

% Check quadrature
assert(abs(exterior_area-4.0*L) < 1e-8)

if gridtype == "quad"
  assert(abs(interior_area-(N-1)*2*L) < 1e-8)
end

% System matrix
Fu = A - S + alpha*S' + P;
condA = condest(Fu)

disp('Solve...')
if alpha == -1
  assert(issymm(Fu))
  [U,flag,relres,it] = pcg(Fu, b, 1e-10, size(Fu,1));
else
  [U,flag,relres,it] = gmres(Fu, b, 10, 1e-10, size(Fu,1));
end

% Display some statistics
relres, it
b_stats = [min(b), max(b), mean(b)]
U_stats = [min(U), max(U), mean(U)]

% Plotting
% Evaluate in midpoint
uh_mid = zeros(G.cells.num, 1);
for e = 1:G.cells.num
  [~, ~, ~, gn] = extractSubgrid(G, e);
  xc = mean(G.nodes.coords(gn, :));
  xhat = disc.transformCoords(xc, e);
  dofs = cell_dofs(e);
  for i = 1:ndof
    uh_mid(e) = uh_mid(e) + dot(U(dofs(i)), psi{i}(xhat));
  end
end
figure
plotCellData(G, uh_mid)
title 'uh\_mid'
uh_mid_stats = [min(uh_mid) max(uh_mid)]

% Solution (approximate) in grid nodes
U_nodal = zeros(G.nodes.num,1);
for e = 1:G.cells.num
  [~, ~, ~, gn] = extractSubgrid(G, e);
  dofs = cell_dofs(e);
  for i = 1:numel(gn)
    x = G.nodes.coords(gn(i),:);
    xhat = disc.transformCoords(x, e);
    for j = 1:ndof
      U_nodal(gn(i)) = U_nodal(gn(i)) + dot(U(dofs(j)), psi{j}(xhat));
    end
  end

end
figure
plotNodeData(G, U_nodal)

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.
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
