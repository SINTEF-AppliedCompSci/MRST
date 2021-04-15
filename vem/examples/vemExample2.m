%% Solving transport problems with VEM
% In this example, we will use the virtual element method (VEM) solver in
% MRST, and show how to set up and use it to solve a transport problem in
% which we inject a high-viscosity fluid in a reservoir with varying
% permeability, which is initially filled with a low-viscosity fluid. To
% emphasize the importance of using a consistent discretization method for
% grids that are not K-orthogonal, we will compare the solution to the
% <matlab:help('incompTPFA') Two-point flux approximation> (TPFA) mehod.

try
   require upr vem incomp vemmech coarsegrid
catch
   mrstModule add upr vem incomp vemmech coarsegrid
end

%% Define geometry
% We will use the UPR module to construct a composite PEBI-grid covering
% the domain $[0, 1000]\times[0, 1000]$, divided into four regions of
% highly varying permeability, with a source placed in the lower left
% corner, and a sink in the upper right corner. See the UPR module on how
% to contruct such grids.
n  = 30;   % Number of cells in each axial direction
l  = 1000; % Physical lenght
% Cell and face constraints
vn = [-1, 1]/sqrt(2);
h = sqrt(2)*l;
f = @(alpha) bsxfun(@plus, bsxfun(@times, vn*min(alpha, 2-alpha), [-h/2; h/2]), [1,1]/sqrt(2)*alpha*h/2);
fc = cellfun(@(alpha) f(alpha), num2cell([0.5, 1, 1.5]), 'UniformOutput', false);
cc = {0.1*l*[1,1], 0.9*l*[1,1]};
% Make grid
G = pebiGrid2D(l/n, [l,l], 'faceConstraints', fc, 'cellConstraints', cc);
% Partition grid into regions
cells = ones(G.cells.num,1);
cells = processPartition(G, cells, find(G.faces.tag));

%%
% Having generated the grid structure, we plot the result.

clf, plotGrid(G, 'facecolor', 'none');
srcCells = find(G.cells.tag);
axis equal off

%%
% The VEM implementation uses grid properties that are not computed by the
% MRST-function <matlab:help('computeGeometry') computeGeometry>, such as
% cell diameters and edge normals. Thus, we compute the geometry using
% <matlab:help('computeVEMGeometry') computeVEMGeometry>. Note that this
% function uses <matlab:help('computeGeometry') computeGeometry>.

G = computeVEMGeometry(G);

%%  Define rock and fluid properties
% Numbering the grid regions from the lower left corner to the upper right
% corner, we set permeability in region $i$ to be
%
% $$ K_i = R(\theta_i)\left(\begin{array}{cc}
%      1000 & 0 \\ 0 & 10
%      \end{array}\right)R(\theta_i)^T, $$
%
% in units of 100 mD, where $R(\theta)$ is a rotation matrix
%
% $$ R(\theta) = \left(\begin{array}{cc}
%      \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta)
%      \end{array}\right). $$
%

K = diag([1000, 10])*milli*darcy();
R = @(t) [cos(t) -sin(t); sin(t) cos(t)];

KK{1} = R(pi/4)*K*R(pi/4)';
KK{2} = K;
KK{3} = R(pi/2)*K*R(pi/2)';
KK{4} = K;

%%
% The structure c consists of four logical vectors identifying the cells
% belonging to each region, and we use this to set the permeability in each
% region. To define a full, symmetric permeability tensor, we specify the
% upper-triangular part. The rock structure is then constructed using
% <matlab:help('makeRock') makeRock>, where we set the porosity to 0.3:

perm = zeros(G.cells.num, 3);
for i = 1:max(cells)
    ix = cells == i;
    perm(ix,:) = repmat(KK{i}([1,2,4]), nnz(ix), 1);
end
    
poro = 0.3;
rock = makeRock(G, perm, poro);

%%
% We will consider a scenario in which we inject a high-viscosity fluid at
% a constant rate in the reservoir, which is initially filled with a fluid
% with lower viscosity. We define a fluid structure with these two fluids
% using <matlab:help('initSimpleFluid') initSimpleFluid>. For simplicity,
% we set both phase relative permeability exponents to 1. Finally, we use
% <matlab:help('initResSol') initResSol> to initialize the state. Here, we
% must indicate that the initial saturation is 0 for the fluid we are
% injecting, and 1 for the fluid which initially occupies the pore volume.

fluid = initSimpleFluid('mu' , [   5,    1]*centi*poise     , ...
                        'rho', [1000,  800]*kilogram/meter^3, ...
                        'n'  , [   1,    1]                      );
                    
[stateVEM, stateTPFA] = deal(initResSol(G, 0, [0,1]));

%% Definig source terms
% We set the injection rate equal one pore volume per 10 years. The source
% and sink terms are constructed using <matlab:help('addSource')
% addSource>. We must specify the saturation in the source term. Since want
% to inject fluid 1, we set this to [1,0]. We must also specify a
% saturation for the sink term, but this does not affect the solution.

Q   = sum(poreVolume(G,rock))/(10*year);
src = addSource([], find(G.cells.tag), [Q, -Q], 'sat', [1,0; 0,1]); 

%% Constructing the linear systems
% Next, we compute the two-point transmissibilties and virtual inner
% products. In this example, we use a first-order VEM.

STPFA = computeTrans(G, rock);
SVEM  = computeVirtualIP(G, rock, 1);

%% Transport
% We now solve the transport problem over a time period equal to the time
% it takes to inject one pore volume. While TPFA is locally conservative by
% construction, VEM is not, which may lead to unphysical saturations when
% the solution is applied to the transport solver. To avoid this, we must
% thus postprocess the VEM solution so that the fluxes are locally
% conservative. This is done using <matlab:help('conserveFlux')
% conserveFlux>, in which we provide any sources, sinks and boundary
% conditions. As for the solvers in MRST, all boundaries without prescribed
% boundary conditions are interpreted as no-flow boundaries.
%
% For each time step, we plot the saturation profiles and the produciton
% rate in the sink cell.

T  = 10*year;
nT = 100;
dT = T/nT;

[productionRateVEM, productionRateTPFA] = deal(zeros(nT,1));

t = linspace(0,T/year,nT);

for i = 1:nT
    
    %   Solve TPFA pressure equation.
    stateTPFA = incompTPFA(stateTPFA, G, STPFA, fluid, 'src', src);
    
    %   Solve TPFA transport equation.
    stateTPFA = implicitTransport(stateTPFA, G, dT, rock, fluid, 'src', src);
    
    %   Plot results.
    subplot(2,2,1)
    plotCellData(G, stateTPFA.s(:,1), 'edgecolor', 'none')
    caxis([0,1]);
    title('TPFA');
    axis equal off
    
    %   Solve VEM pressure equation.
    stateVEM = incompVEM(stateVEM, G, SVEM, fluid, 'src', src, 'cellPressure', true);
    
    %   Postprocess solution.
    stateVEM = conserveFlux(stateVEM, G, rock, 'src', src);
    
    %   Solve VEM transport equation.
    stateVEM = implicitTransport(stateVEM, G, dT, rock, fluid, 'src', src);
    
    %   Plot the results.
    subplot(2,2,2)
    plotCellData(G, stateVEM.s(:,1), 'edgecolor', 'none')
    caxis([0,1]);
    title('VEM');
    axis equal off
    
    %   Plot production in producer.
    productionRateVEM(i)  = day*Q*stateTPFA.s(src.cell(2),2);
    productionRateTPFA(i) = day*Q*stateVEM.s(src.cell(2),2);
    
    subplot(2,2,3:4)
    h = plot(t(1:i), productionRateVEM( 1:i), ...
             t(1:i), productionRateTPFA(1:i), 'lineWidth', 2);
    axis([0 T/year 0 1.1*day*Q])
    xlabel('Time [years]');
    ylabel('Production rate [m^3/day]');
    legend('TPFA', 'VEM')
    drawnow
    
end

%%
% We see that there are significant differences in the two saturation
% profiles and production curves. This is due to the fact that TPFA fails
% to capture the effect of the rotated permeability fields in regions one
% and four.

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
