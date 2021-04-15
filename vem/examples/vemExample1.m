%% Basic VEM tutorial
% In this example, we will introduce the virtual element method (VEM)
% solver in MRST, and show how to set up and use it to solve the
% single-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,
%    \qquad x \in \Omega,$$
%
% for a simple pressure drop problem.
%
% The virtual element method constitutes a unified framework for
% higher-order methods on general polygonal and polyhedral grids. We
% combine the single-phase pressure equaiton into an ellpitic equation in
% the pressure only;
%
% $$ - \nabla \cdot K\nabla p = q.$$
%
% Multiplying this by a function $v \in H_0^1(\Omega)$ and integrating over
% $\Omega$ yields the weak formulation: Find $p \in H^1_0(\Omega)$ such
% that
%
% $$ a(p, v) = \int_{\Omega}\nabla p \cdot K \nabla v d x
%   = \int_\Omega q v d x, \qquad \forall v \in H^1_0(\Omega).$$
%
% The bilinear form $a$ can be split into a sum of local bilinear forms as
%
% $$ a(u,v) = \sum\nolimits_i
%     \int_{\Omega_i} \nabla u \cdot K \nabla v d x
%   = \sum\nolimits_i a^i(u,v),$$
%
% where $\Omega_i$ are the cells of the grid. For a $k$-th order VEM, we
% construct approximate local bilinear forms $a_h^i$ which are exact
% whevener one of $u,v$ is a polynomial of degree $k$ or less. When none of
% $u,v$ are a such a polynomial, we only approximate $a_h^i$ to the right
% order. To obtain this, we apply the following Pythagoras identity:
%
% $$ a^i(u,v) = a^i(\Pi u, \Pi v) + a^i(u-\Pi u, v- \Pi v),$$
%
% where $\Pi$ is an $a^i$-orthogonal projection onto the space of
% polynomials of degree $k$ or less. The approximate bilinear form is then
%
% $$ a_h^i(u,v) = a^i(\Pi u, \Pi v) + s^i(u-\Pi u, v- \Pi v),$$
%
% where $s^i$ is a stability term, which we are free to choose so long as
% the resulting bilinear form is positive definite and stable.  The flux
% field can easily be reconstructed from the pressure solution. For details
% on constructing the VEM function spaces, and implementing the method, see
% for example <http://hdl.handle.net/11250/2405996 Klemetsdal 2016>.

try
   require upr vem incomp vemmech
catch
   mrstModule add upr vem incomp vemmech
end

%% Define geometry
% Since VEM is consistent for general polyhedral cells, we will use the UPR
% module to construct an arbitrary PEBI-grid. We consider a domain covering
% $[0, 500]\times[0, 100]\times[0, 100]$. To generate the grid, will use
% the function <matlab:help('mirroredPebi') mirroredPebi>, in which we
% define the boundary of the domain, and the seed-points for the gridcells.

gridLimits = [500, 100, 100];

boundary  = [0             0             0             ; ...
             gridLimits(1) 0             0             ; ...
             gridLimits(1) gridLimits(2) 0             ; ...
             0             gridLimits(2) 0             ; ...
             0             0             gridLimits(3) ; ...
             gridLimits(1) 0             gridLimits(3) ; ...
             gridLimits(1) gridLimits(2) gridLimits(3) ; ...
             0             gridLimits(2) gridLimits(3)];

points = bsxfun(@times, rand(200, 3), gridLimits);

G = mirroredPebi3D(points, boundary);

%%
% Having generated the grid structure, we plot the result.

clf, plotGrid(G); view(3); axis equal
light('Position',[-50 -50 -100], 'Style', 'infinite')

%%
% The VEM implementation uses grid properties that are not computed by the
% MRST-function <matlab:help('computeGeometry') computeGeometry>, such as
% cell diameters and edge normals. Thus, we compute the geometry using
% <matlab:help('computeVEMGeometry') computeVEMGeometry>. Note that this
% function uses <matlab:help('computeGeometry') computeGeometry> and the
% function <matlab:help('createAugmentedGrid') createAugmentedGrid> from
% the VEM module for geomechanics.

G = computeVEMGeometry(G);

%%  Define rock and fluid properties
% We set the permeability to be homogeneous and anisotropic,
%
% $$ K = \left(\begin{array}{ccc}
%      1000 & 0  & 0 \\ 0 & 100 & 0 \\ 0 & 0 & 10 \end{array}\right). $$
%
% The rock structure is constructed usin <matlab:help('makeRock')
% makeRock>, where we let the porosity be 0.3 in all cells. The fluid is
% constructed using <matlab:help('initSinglefluid') initSingleFluid>.
% Finally, we use <matlab:help('initResSol') initResSol> to initialize the
% state.

perm = [1000, 100, 10].* milli*darcy();
poro = 0.3;
rock = makeRock(G, perm, poro);

fluid = initSingleFluid('mu' , 100*centi*poise     , ...
                        'rho', 800*kilogram/meter^3);

state = initResSol(G, 0);

%% Definig boundary conditions
% Since the PEBI-grid might be subject to roundoff-errors on the
% boundaries, we identify specific boundary faces by finding the ones that
% are within a given tolerance away from where we expect the boundary to
% be.

tol = 1e-6;
                    
bFaces = boundaryFaces(G);
west   = abs(G.faces.centroids(bFaces,1)                ) < tol;
east   = abs(G.faces.centroids(bFaces,1) - gridLimits(1)) < tol;

%%
% We define the western and easter boundaries to be Dirichlet boundaries
% with a prescribed pressure. Boundaries with no prescribed boundary
% condition are no-flow boundaries by default.

bc = addBC([], bFaces(west), 'pressure', 500*barsa);
bc = addBC(bc, bFaces(east), 'pressure', 0        );

%% Constructing the linear systems
% First, we compute the local inner products $a^i$. The implementation
% supports first and second order, and we will compute both. This is done
% using <matlab:help('computeVirtualIP') computeVirtualIP>, where we must
% provide the grid and rock structures, and method order.

S1 = computeVirtualIP(G, rock, 1); disp(S1);
S2 = computeVirtualIP(G, rock, 2);

%%
% The resulting solution structure has the following fields:
%
% # A: Block diagonal matrx with the local bilinear forms on the
%   diagonal.
% # ip: Choice of stability term for the inner products. Possible choices
% are 'ip_simple', 'ip_qfamily', 'ip_fem', and 'ip_fd'.
% # dofVec: Map from local to global degrees of freedom.
% # PiNstar and PiNFstar: Matrix representations of the projection
% operator $\Pi$ on the cells and faces of the grid.
% # faceCoords: Local 2D coordinate systems on each face.
% # order: Method order.
% # T, transType: The implementation also supports flux reconstruction
% using the TPFA or MPFA schemes. This can be obtained by setting the
% 'trans' option to 'tpfa' or 'mpfa'. In this case, T and transType holds
% the half-face transmissibilities, and the corresponding scheme.
%
% Next, we compute the solution using both first-order and second-order
% VEM. The first-order VEM uses the pressure at the nodes as degrees of
% freedom. However, it is possible to reconstruct the cell pressures of the
% approximated solution exactly, by uising the node pressures. This can be
% done by setting the optional input parameter 'cellPressure' to true in
% <matlab:help('incompVEM') incompVEM>. Moreover, we will invesitgate the
% linear system of the two methods later, so we set 'matrixOutput' to true.

state1 = incompVEM(state, G, S1, fluid, 'bc', bc, ...
                            'cellPressure', true, 'matrixOutput', true);
state2 = incompVEM(state, G, S2, fluid, 'bc', bc, 'matrixOutput', true);

%%
% We then plot the results.

clf

subplot(2,1,1)
plotCellData(G, state1.pressure); axis equal off; colorbar; view(3);

subplot(2,1,2)
plotCellData(G, state2.pressure); axis equal off; colorbar; view(3);

%% Degrees of freedom
% Finally, we investigate the linear systems obtained from the two methods.
% We start by identifying the cell closest to the middle of the domain, its
% node coordinates, and face and cell centroids.

distances = sum(bsxfun(@minus, G.cells.centroids, gridLimits/2).^2,2);
midCell   = find(distances == min(distances));
midCell   = midCell(1);

n   = G.cells.nodes(G.cells.nodePos(midCell):G.cells.nodePos(midCell+1)-1);
xn  = G.nodes.coords(n,:);
f   = G.cells.faces(G.cells.facePos(midCell):G.cells.facePos(midCell+1)-1);
xf  = G.faces.centroids(f,:);
xc  = G.cells.centroids(midCell,:);

%%
% The first-order VEM uses the pressure at the nodes as degrees of freedom.
% The second-order VEM uses the pressure at the nodes and edge centroids,
% along with the average face pressure and average cell pressure. The
% latter yields a much denser system. We show this by indicating the
% degrees of freedom of the middle cell, along with a sparsity plot of the
% linear systems.

clf
subplot(2,2,1)
plotGrid(G, midCell, 'faceAlpha', 0.2); axis equal off;
hold on
plot3(xn(:,1), xn(:,2), xn(:,3), 'sq');

subplot(2,2,2)
spy(state1.A)

subplot(2,2,3)
plotGrid(G, midCell, 'faceAlpha', 0.2); axis equal off;
hold on
plot3(xn(:,1), xn(:,2), xn(:,3), 'sq')
plot3(xf(:,1), xf(:,2), xf(:,3), 'd')
plot3(xc(:,1), xc(:,2), xc(:,3), 'o')

subplot(2,2,4)
spy(state2.A)

%%
% The construciton of these linear systems can be quite expensive in terms
% of computational effort, especially for reservoir models with many cells.
% In order to speed up the construction, one can use C-accelerated MEX
% functions in <matlab:help('computeVirtualIP') computeVirtualIP>. This is
% done by setting 'invertBlocks' to 'MEX'. Note that this relies on being
% able to build the required MEX functions:

tic
fprintf('Computing innerproducts using mldivide \t ... ');
computeVirtualIP(G, rock, 2);
toc

tic
fprintf('Computing innerproducts using MEX \t ... ');
computeVirtualIP(G, rock, 2, 'invertBlocks', 'MEX');
toc

%%
% For larger resrvoir models, using C-accelerated MEX functions is usually
% several orders of magnitude faster, and should always be used if
% possible.

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
