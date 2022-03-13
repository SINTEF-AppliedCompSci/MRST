%% Grid orientation example
% This example is a slightly revised version of the grid-orientation
% example from Section 10.4.2 of the MRST textbook. The fluid system is
% incompressible, and to be able to simulate the system using AD-OO solvers
% we add a small amount of rock compressibility. The pressure solver will
% complain about a badly scaled matrix but will nonetheless give an
% acceptable pressure solution.

%% Add modules
mrstModule add dg ad-core ad-props ad-blackoil sequential

%%
gravity reset off
pvi      = 0.275;
cartDim  = [32 32];
nstep    = 64;
T        = year;
fluid    = initSimpleADIFluid('phases', 'WO' , 'n', [2,3], ... 
                              'mu', [1,10],  'rho', [1,1], 'cr', 1e-8/barsa);

%% Original setup
domain = [1000 1000];
G      = computeGeometry(cartGrid(cartDim,domain));
G      = computeCellDimensions(G);
rock   = makeRock(G, 450*milli*darcy, 0.2);
rate   = pvi*sum(poreVolume(G,rock))/year;
W      = addWell([],G, rock, 1, 'Type', 'rate', ...
          'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W      = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
          'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
model  = GenericBlackOilModel(G, rock, fluid, 'gas', false);
sched  = simpleSchedule(rampupTimesteps(T,T/nstep), 'W', W);

%% Rotated quater-five spot
theta  = pi/4;
R      = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Gr     = cartGrid(round(G.cartDims.*sqrt(2)), domain);
Gr.nodes.coords = sqrt(2)*(R*(Gr.nodes.coords'))';
Gr     = computeGeometry(Gr);
Gr      = computeCellDimensions(Gr);
rockr  = makeRock(Gr, 450*milli*darcy, 0.2);
n      = Gr.cartDims(1);
Wr     = addWell([], Gr, rockr, 1, 'Type', 'rate', ...
          'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wr     = addWell(Wr, Gr, rockr, n*n, 'Type', 'rate', ...
          'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
Wr     = addWell(Wr, Gr, rockr, n, 'Type', 'rate', ...
          'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
Wr     = addWell(Wr, Gr, rockr, n*(n-1)+1, 'Type', 'rate', ...
          'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
modelr = GenericBlackOilModel(Gr, rockr, fluid, 'gas', false);
schedr = simpleSchedule(rampupTimesteps(T,T/nstep), 'W', Wr);

%% Finite volume simulations
pmodel  = PressureModel (model);
tmodel  = TransportModel(model);
dmodel  = TransportModelDG(model, 'degree', [1,1]);
pmodelr = PressureModel (modelr);
tmodelr = TransportModel(modelr);
dmodelr  = TransportModelDG(modelr, 'degree', [1,1]);
[~, x ] = simulateScheduleAD(initResSol(G,  100*barsa, [0,1]), ...
    SequentialPressureTransportModel(pmodel, tmodel, 'parentModel', model), sched);
[~, xr] = simulateScheduleAD(initResSol(Gr, 100*barsa, [0,1]),...
    SequentialPressureTransportModel(pmodelr, tmodelr, 'parentModel', modelr), schedr);

%%
[~, d ] = simulateScheduleAD(initResSol(G,  100*barsa, [0,1]), ...
    SequentialPressureTransportModel(pmodel, dmodel, 'parentModel', model), sched);
[~, dr] = simulateScheduleAD(initResSol(Gr, 100*barsa, [0,1]),...
    SequentialPressureTransportModel(pmodelr, dmodelr, 'parentModel', modelr), schedr);


%% Plot the results
N        = 20;
cval     = [.5*movsum(linspace(0,1,N+1),2) 1];
plotData = @(G,x)...
    contourf(reshape(G.cells.centroids(:,1),G.cartDims), ...
             reshape(G.cells.centroids(:,2),G.cartDims), ...
             reshape(x,G.cartDims), cval,'EdgeColor','none');
Gp = cartGrid(Gr.cartDims.*[1 2], domain.*[1 2]);
Gp.nodes.coords(:,2) = Gp.nodes.coords(:,2) - domain(2);
Gp.nodes.coords = sqrt(2)*(R*(Gp.nodes.coords'))';
Gp = computeGeometry(Gp);

colormap(flipud([.7*winter(128).^2+.3; 1 1 1]));

for n=1:2
    
    subplot(1,2,n)
    if n==1
        s  = x{end-1} .s(:,1);
        sr = xr{end-1}.s(:,1);
    else
        s  = d{end-1} .s(:,1);
        sr = dr{end-1}.s(:,1);
    end
    
    plotData(G,s)
    hold on
    sr = reshape(sr, Gr.cartDims); sr = sr(:,[end:-1:1 1:end],:);
    contour(reshape(Gp.cells.centroids(:,1), Gp.cartDims),...
        reshape(Gp.cells.centroids(:,2), Gp.cartDims), ...
        reshape(sr(:),Gp.cartDims), cval(2:end-1),'k');
    hold off
    axis equal tight; axis([0 domain(1) 0 domain(2)]);
end

%% Copyright Notice
%
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
