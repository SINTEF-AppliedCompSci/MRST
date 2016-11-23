%% Solving transport problems with VEM
% In this example, we will use the virtual element method (VEM) solver in
% MRST, and show how to set up and use it to solve the single-phase
% pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\frac{K}{\mu}\nabla p,$$
%
% for a problem with one source and one sink. We will then use the
% resulting solution to solve the associated transport problem. To
% emphasize the importance of using a consistent discretization method for
% grids that are not K-orthogonal, we will compare the solution to the
% Two-point flux approximation (TPFA) mehod.

try
   require upr vem incomp
catch
   mrstModule add upr vem incomp
end

%% Define geometry
% We will use the UPR module to construct fully unstructured
% PEBI-grid covering the domain $[0, 1000]\times[0, 1000]$. To generate the
% grid, will use the function <matlab:help('pebiGrid') pebiGrid>, in which
% we define the average gridcell size, the physical grid sizes in the axial
% directions, and the well coordinates.

gridLimits = [1, 1];
wellCoordinates = {[0.1, 0.1].*gridLimits; [0.9, 0.9].*gridLimits};
n = 20;

G = pebiGrid(gridLimits(1)/n, gridLimits, 'wellLines', wellCoordinates);

G.nodes.coords = G.nodes.coords*1000;
G = computeVEMGeometry(G);

%%
% Having generated the grid structure, we plot the result. The source and
% sink cells are indicated in blue and red, respecitvely.

clf
plotGrid(G);
srcCells = find(G.cells.tag);
plotGrid(G, srcCells(1), 'facecolor', 'b');
plotGrid(G, srcCells(2), 'facecolor', 'r');
axis equal off

%%
% The VEM implementation uses grid properties that are not computed by the
% MRST-function <matlab:help('computeGeometry') computeGeometry>, such as
% cell diameters and edge normals. Thus, we compute the geometry using
% <matlab:help('computeVEMGeometry') computeVEMGeometry>. Note that this
% function uses <matlab:help('computeGeometry') computeGeometry>.

G = computeVEMGeometry(G);

%%  Define rock and fluid properties
% We set the permeability to be homogeneous and anisotropic
%
% $$ K = R(\pi/4)\left(\begin{array}{cc}
%      1000 & 0 \\ 0 & 10
%      \end{array}\right)R(\pi/4)^T, $$
%
% where $R(\theta)$ is a rotation matrix
%
% $$ R(\theta) = \left(\begin{array}{cc}
%      \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta)
%      \end{array}\right). $$
%
% We use <matlab:help('makeRock') makeRock>, to define this permeability
% tensor, and let the porosity be 0.3 in all cells. We will consider a
% scenario in which we inject water at a constant rate in the reservoir,
% which is initially filled with oil. We define a fluid structure with
% these two fluids using  <matlab:help('initSimpleFluid') initSimpleFluid>.
% For simplicity, we set both Corey coefficients to 1. Finally, we use
% <matlab:help('initResSol') initResSol> to initialize the state. Here, we
% must indicate that the initial saturation is 1 for oil, and 0 for water.

perm = diag([1000, 10])*milli*darcy();
R = @(t) [cos(t) -sin(t); sin(t) cos(t)];
t = pi/4;
perm = R(t)*perm*R(t)';

poro = 0.3;
rock = makeRock(G, perm([1,2,4]), poro);

fluid = initSimpleFluid('mu' , [   1,  100]*centi*poise     , ...
                        'rho', [1014,  800]*kilogram/meter^3, ...
                        'n'  , [   1,    1]                      );
                    
[stateVEM, stateTPFA] = deal(initResSol(G, 0, [0,1]));

%% Definig source terms
% We assume that our domain extends 100 m in the $z$-direction, and set the
% injection rate equal one pore volume per 10 years. The source and sink
% terms are constructed using <matlab:help('addSource') addSource>. For the
% source term, we must specify the saturation of the injected fluid. Since
% want to inject water, we set this to [1,0].

Q = 100*sum(poreVolume(G,rock))/(10*year);
src = addSource([], find(G.cells.tag), [Q, -Q], 'sat', [1,0; 0,1]); 

%% Constructing the linear systems
% We compute the two-point transmissibilties and virtual inner products. In
% this example, we use a first-order VEM.

STPFA = computeTrans(G, rock);
SVEM  = computeVirtualIP(G, rock, 1);

%% Transport
% We now solve the transport problem over a time period equal to the time
% it takes to inject one pore volume. While TPFA is locally conservative by
% construction, VEM is not, which may lead to unphysical saturations when
% the solution is applied to the transport solver. To avoid this, we must
% thus postprocess the VEM solution. This is done using
% <matlab:help('conserveFlux') conserveFlux>, in which we provide any
% sources, sinks and boundary conditions. as for the MRST-solver, all
% boundaries without prescribed boundary conditions are interpreted as
% no-flow boundaries.
%
% For each time step, we plot the saturation profiles and the produciton
% rate in the sink cell.

T = 10*year;
nT = 50;
dT = T/nT;

[productionRateVEM, productionRateTPFA] = deal(zeros(nT,1));

t = linspace(0,T/year,nT);
r = cell(nT,1);

for i = 1:nT
    
    stateTPFA = incompTPFA(stateTPFA, G, STPFA, fluid, 'src', src);
    stateTPFA = implicitTransport(stateTPFA, G, dT, rock, fluid, 'src', src);
    
    subplot(2,2,1)
    plotCellData(G, stateTPFA.s(:,1), 'edgecolor', 'none')
    caxis([0,1]);
    title('TPFA');
    axis equal off
    
    stateVEM = incompVEM(stateVEM, G, SVEM, fluid, 'src', src);
    [stateVEM, r{i}] = conserveFlux(stateVEM, G, rock, 'src', src);
    stateVEM = implicitTransport(stateVEM, G, dT, rock, fluid, 'src', src);
    
    subplot(2,2,2)
    plotCellData(G, stateVEM.s(:,1), 'edgecolor', 'none')
    caxis([0,1]);
    title('VEM');
    axis equal off
    
    productionRateVEM(i)  = day*Q*stateTPFA.s(src.cell(2),2);
    productionRateTPFA(i) = day*Q*stateVEM.s(src.cell(2),2);
    
    subplot(2,2,3:4)
    plot(t(1:i), productionRateVEM(1:i), t(1:i), productionRateTPFA(1:i))
    axis([0 T/year 0 day*Q])
    xlabel('Time [years]');
    ylabel('Production rate [m^3/day]');
    legend('TPFA', 'VEM')
    pause(0.0001)
    
end