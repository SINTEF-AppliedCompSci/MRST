%% Introduction to Flow Diagnostics
% In this tutorial, we will introduce you to the basic ideas of flow
% diagnostics. The basic quantitites in flow diagnostics are time-of-flight
% and numerical tracer partitions, which in turn can be used to compute
% volumetric connections, flux allocation between wells, and measures of
% dynamic heterogeneity, or just be used to indicate how displacement
% fronts will propagate in the reservoir, given a steady flow field.
%
% Flow diagnostics has often been associated with streamline methods.
% Herein, however, we will compute the same kind of quantities and measures
% using the first-order, finite-volume discretization which is used in MRST
% and in most commercial reservoir simulators. This means that whereas the
% computed flow diagnostics may not necessarily be an accurate
% representation of the corresponding continuous problem, it will represent
% the volumetric connections and dynamic heterogeneity, *as seen by a
% multiphase simulator*.
%
%
% Suggested reading:
% 
% # O. Møyner, S. Krogstad, and K.-A. Lie. The application of flow
% diagnostics for reservoir management. SPE J., Vol. 20, No. 2, pp.
% 306-323, 2015. DOI: 10.2118/171557-PA
% #  M. Shahvali, B. Mallison, K. Wei, and H. Gross. An alternative to
% streamlines for flow diagnostics on structured and unstructured grids.
% SPE J., 17(3):768{778, 2012. DOI: 10.2118/146446-PA
% # Shook, G.M. and Mitchell, K.M. 2009. A robust measure of heterogeneity
% for ranking earth models: The F PHI curve and dynamic Lorenz coefficient.
% SPE Annual Technical Conference and Exhibition, New Orleans, 4-7 October.
% DOI: 10.2118/124625-MS.

mrstModule add diagnostics incomp streamlines;

% NB: In this example, we use the |parula| colormap. If this is not yet
% available in your MATLAB, you can e.g., replace all occurences of
% |parula| with |jet|
pargs = {'EdgeColor','none','FaceAlpha',.6};

%% Set up and solve flow problem
% To illustrate the various concepts, we use a rectangular reservoir with
% five wells, two injectors an five producers

% Grid
[nx,ny] = deal(64);
G = cartGrid([nx,ny,1],[500,250,10]);
G = computeGeometry(G);

% Petrophysical data
p = gaussianField(G.cartDims(1:2), [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1.5e-5)^2./(0.81*72*(1-p).^2);

rock = makeRock(G, K(:), p(:));

hT  = computeTrans(G, rock);
pv  = sum(poreVolume(G,rock));

% Fluid model
gravity reset off
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

% Wells
n = 12;
W = addWell([],  G, rock, nx*n+n+1, ...
    'Type', 'rate', 'Comp_i', 1, 'name', 'I1', 'Val', pv/2);
W = addWell(W, G, rock, nx*n+n+1+nx-2*n, ...
    'Type','rate',  'Comp_i', 1, 'name', 'I2', 'Val', pv/2);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx-.3*nx), ...
    'Type','rate',  'Comp_i', 0, 'name', 'P1', 'Val', -pv/4);
W = addWell(W, G, rock, G.cells.num-(n-.5)*nx, ...
    'Type','rate',  'Comp_i', 0, 'name', 'P2', 'Val', -pv/2);
W = addWell(W, G, rock, round(G.cells.num-(n/2-.5)*nx+.3*nx), ...
    'Type','rate',  'Comp_i', 0, 'name', 'P3', 'Val', -pv/4);

% Initial reservoir state
state = initState(G, W, 0.0, 1.0);

%% Compute basic quantities
% To compute the basic quantities, we need a flow field, which we obtain by
% solving a single-phase flow problem. Using this flow field, we can
% compute time-of-flight and numerical tracer partitions. For comparison,
% we will also tracer streamlines in the same flow field.

state = incompTPFA(state, G, hT, fluid, 'wells', W);
D = computeTOFandTracer(state, G, rock, 'wells', W);

% Trace streamlines
seed = (nx*ny/2 + (1:nx)).';
Sf = pollock(G, state, seed, 'substeps', 1);
Sb = pollock(G, state, seed, 'substeps', 1, 'reverse', true);

%% Forward time-of-flight
% The time-of-flight gives the time it takes a neutral particle to travel
% from the nearest inflow boundary (e.g., an injection well) and to a given
% point in the reservoir. Here, we do not compute this quantity pointwise,
% but rather in an volume-averaged manner over each grid cell
clf
hf=streamline(Sf);
hb=streamline(Sb);
plotGrid(G,vertcat(W.cells),'FaceColor','none','EdgeColor','r','LineWidth',1.5);
plotWell(G,W,'FontSize',20); 
set([hf; hb],'Color','w','LineWidth',1.5);
hd=plotCellData(G,D.tof(:,1),pargs{:});
colormap(parula(32));
axis equal off;

%% Backward time-of-flight
% Backward time-of-flight is defined analogously, as the name suggests, by
% tracing particles backward along the flow field from the nearest outlet
% (e.g., a production well). For a given point, this quantity thus
% represents the time it will take a neutral particle to reach the outlet.
delete(hd);
hd=plotCellData(G,D.tof(:,2),pargs{:});

%% Residence time
% The residence time is defined as the sum of the forward and backward
% time-of-flight and gives the total travel time from inlet to outlet. This
% quantity can be used to distinguish high-flow zones, which have small
% residence time, from low-flow and stagnation zones, which have high
% residence times.
delete(hd);
hd=plotCellData(G,sum(D.tof,2),pargs{:});

%% Tracer from I1
% Numerical tracer partitions can be thought of as an imaginary paiting
% experiment, in which we inject a massless, nondiffusive ink of a unique
% color at each fluid source. At time infinity, each point in the reservoir
% must then have been colored. If a point receives more than one tracer,
% the tracer concentrations represent the fractional flow of past each
% points. The result is a volumetric partition of unity for the reservoir.
% To illustrate, let us look at the tracer concentration of injector I1,
% which is one in all parts of the reservoir except near the flow divide,
% where cells are affected by flow from both injectors.
delete(hd);
t = D.itracer(:,1);
hd = plotCellData(G, t, t>1e-6, pargs{:});

%% Flooded regions
% Using this injection tracer, we can define the flooded regions. Each
% inejctor will influence all cells having a positive tracer concentration.
% In many cases, however, it is interesting to determine the injector that
% influences the most. In MRST, this is done by a simple majority vote over
% all tracers. In lack of a better name, this will be referred to as 'flood
% regions'
delete(hd);
hd = plotCellData(G,D.ipart,pargs{:});

%% Tracer from P2
% By reversing the flow field, we can compute a similar partition of unity
% associated with the producers, here illustrated for P2, In the figure, we
% see that P2 has a region in the middle of the domain that it drains all
% by itself, whereas the regions to the east and west are also drained by
% P3 and P1
delete(hd);
t = D.ptracer(:,2);
hd = plotCellData(G, t, t>1e-6, pargs{:});

%% Drainage regions
% Like we associated flood regions with injectors, we can associate
% drainage regions with producers, using a simple majority vote over all
% tracer concentrations
delete(hd);
hd = plotCellData(G,D.ppart,pargs{:});

%% Well regions: I1<->P1, I2<->P3
% By combining the drainage and flood regions, we can get regions that are
% affected by a single injector-producer pair
delete(hd);
hd = plotCellData(G,D.ipart, D.ppart~=2, pargs{:});

%% Compute time-of-flights inside each well region
% If a cell is affected by more than one injector (or producer), the volume
% averaged time-of-flight can be very inaccurate as a pointwise measure if
% the travel times from the two injectors are very different. To amend
% this, we can compute time of flight for each well
T = computeTimeOfFlight(state, G, rock, 'wells', W, ...
    'tracer',{W(D.inj).cells},'computeWellTOFs', true);

%% F-Phi diagram
% To define a measure of dynamic heterogeneity, we can think of the
% reservoir as a bundle of non-coummunicating volumetric flow paths
% (streamtubes) that each has a volume, a flow rate, and a residence time.
% For a given time, the storage capacity Phi, is the fraction of flow paths
% in which fluids have reached the outlet, whereas F represent the
% corresponding fractional flow. Both are monotone functions of residence
% time, and by plotting F versus Phi, we can get a visual picture of the
% dynamic heterogeneity of the problem. In a completely homogeneous
% displacement, all flowpaths will break through at the same time and hence
% F(Phi) is a straight line from (0,0) to (1,1). In a heterogeneous
% displacement, F(Phi) will be a concave function in which the steep
% initial slope corresponds to high-flow regions giving early breakthrough
% and, whereas the flat trailing tail corresponds to low-flow and stagnant
% regions that would only break through after very long time
[F,Phi] = computeFandPhi(poreVolume(G,rock), D.tof);
clf, plot(Phi,F,'.',[0 1],[0 1],'--');

%% Lorenz coefficient
% The further the F(Phi) curve is from a straight line, the larger the
% difference between fast and slow flow paths. The Lorenz coefficient
% measures the dynamic heterogeneity as two times the area between F(Phi)
% and the straight line F=Phi.
computeLorenz(F,Phi)

%% Sweep effciency diagram
% We can also define measures of sweep efficiency that tell how effective
% injected fluids are being used. The volumetric sweep efficiency Ev is
% defined as the ratio of the volume that has been contacted by the
% displacing fluid at time t and the volume contacted at infinite time.
% This quantity is usually related to dimensionless time td=dPhi/dF
[Ev,tD] = computeSweep(F,Phi);
clf, plot(tD,Ev,'.');

%% Copyright notice

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
