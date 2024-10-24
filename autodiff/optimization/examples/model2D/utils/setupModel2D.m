% setupModel:  construct model for example analyseModel2D

% set up grid and rock-properties for simple 2D-model
% nxyz = [40, 40, 1];
% Dxyz = [400, 400, 10];
% 
% G = cartGrid(nxyz, Dxyz);
% G = computeGeometry(G);
% 
% rock = getSPE10rock(1:40, 101:140, 2);
clear all;
cartDims = [  60,  60,  1];
physDims = [600, 600, 10];
G  = cartGrid(cartDims, physDims);
G  = computeGeometry(G);
[T,layer] = deal(10*year, [23, 75]);
n = 2;
 rock = getSPE10rock(1:cartDims(1),1:cartDims(2),layer(n));
rock.poro = max(rock.poro, 0.14);
% fluid
pRef = 200*barsa;

fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 1e-5/barsa;
p_ref = 200*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

W = [];
% Injectors (lower-left and upper-right)
ci(1) = 1;
ci(2) = G.cells.num;
for k  = 1:2
    W = addWell(W, G, rock, ci(k), 'Type' , 'rate', ...
                                   'Val'  , 300*meter^3/day, ...
                                   'Name' , sprintf('I%d', k), ...
                                   'comp_i', [1 0], ...
                                   'Sign' , 1);
end
% Producers (upper-left and -right)
cp(1) = G.cartDims(1);
cp(2) = 1 + (G.cartDims(2)-1)*G.cartDims(1);
for k  = 1:2
    W = addWell(W, G, rock, cp(k), 'Type', 'bhp', ...
                                   'Val' , 150*barsa, ...
                                   'Name', sprintf('P%d', k), ...
                                   'comp_i', [0 1], ...
                                   'Sign', -1);
end


% Set up 4 control-steps each 150 days
ts = { [0.1 0.25 0.5 1 3 3 3 5 5 5 5 5 5 5 5 5]'*day, ...
                   repmat(70/10, 10, 1)*day, ...
                   repmat(80/10, 10, 1)*day};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:3)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

gravity on


% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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

