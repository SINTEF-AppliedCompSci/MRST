%% Some pressure solvers may include the effect of regions, by the hack
%% "call fluid.relperm with as many saturations as there are cells in G".
nx=100;
G = cartGrid([nx,1], [100,10]);
G = computeGeometry(G);

rock.perm  = ones(G.cells.num, 1);%logNormLayers(G.cartDims, 'a', 3);
rock.poro  = ones(G.cells.num, 1);
rock.poro(G.cells.centroids(:,1)<22 & G.cells.centroids(:,1)>18) = 1e-3;
%% Make special fluid:
satnum = ones(G.cartDims);
%satnum(1:end/2, end/2+1:end)     = 2;
%satnum(end/2+1:end, 1:end/2)     = 3;
%satnum(end/2+1:end, end/2+1:end) = 4;
%satnum = satnum(:);

fluid      = initCoreyFluidROC('mu',  [1,1], ...
                               'rho', [0,0], ...
                               'sr',  [0,0;  ...
                                       0,0;  ...
                                       0,0;  ...
                                       0,0], ...
                               'n',   [1,1;  ...
                                       2,2;  ...
                                       2,3;  ...
                                       2,4], ...
                               'kwm', [1,1;  ...
                                       1,1;  ...
                                       1,1;  ...
                                       1,1], ...
                               'reg', satnum);

%% Compute pressure...
% Compute transport
cfl_vec=[1000,100,10,1,0.5,0.2,0.1];
T          = computeTrans(G, rock);
ss=[];
for kk=1:numel(cfl_vec)
disp('Computing transmissibilities...')

%S          = computeMimeticIP(G, rock);
state      = initState(G, [], 0);
src        = addSource([], [1, G.cells.num], [1, -1], 'sat', [1,1]);
disp('Computing pressure...')
state      = incompTPFA(state, G, T, fluid, 'src', src);
%state      = solveIncompFlow(state, G, S, fluid, 'src', src);



nn=10*cfl_vec(kk);
dt = repmat(10/nn, [nn, 1]);
disp('Computing transport...')
cla,
xc=G.cells.centroids(:,1);

for i=1:70,
   state = implicitupwind(state, G, dt, rock, fluid, 'src', src);
   %plotCellData(G, state.s(:,1)); axis equal tight
   plot(xc,state.s)
   pause(0.01)
   %surf(reshape(state.s(:,1), G.cartDims));shading flat;
   %set(gca, 'dataaspectr', [1,1,0.02]);view(130, 10)
   %drawnow;
end
ss=[ss,state.s(:,1)];
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.
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
