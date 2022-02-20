%% Compare Time-of-Flight Computed by Steamlines and by Finite-Volumes
% MRST offers two ways of computing time-of-flight for Cartesian grids:
% 
% * using Pollock tracing from the |streamlines| module
% * using finite-volume methods from the |diagnostics| module
% 
% Pollock's method gives pointwise values of high accuracy for Cartesian
% grids, but this method has not yet been implemented for general
% unstructured grids. The finite-volume method computes time-of-flight in a
% volume-averaged sense, and will hence not have high pointwise accuracy,
% but can readily be applied to any type of valid MRST grid.
mrstModule add mimetic incomp diagnostics streamlines

%% Setup model and solve flow problem
% We consider a 60x60 grid with a quarter five-spot well pattern. The
% petrophysical quantities can either be uniform, or sampled from a
% Gaussian distribution. Pollock's method is not yet implemented for
% non-unit porosity values, but this minor deficiency can be circumvented
% by appropriate postprocessing and scaling of cell increments

G = cartGrid([60,60]);
G = computeGeometry(G);

%perm = 10*milli*darcy;
perm = convertFrom(logNormLayers([G.cartDims, 1], 1), milli*darcy);
poro = 1;
rock = makeRock(G, perm, poro);

fluid = initSimpleFluid('mu',  [   1,    1]*centi*poise, ...
                        'rho', [1000, 1000]*kilogram/meter^3, ...
                        'n',   [   2,    2]);

IP = computeMimeticIP(G, rock);
src = addSource([], [1, G.cells.num], [1.1, -1.1], 'sat', [1.0 0.0;1.0,0.0]);

x = initResSol(G, 0, 0);
x = incompMimetic(x, G, IP, fluid, 'src', src);

%% Trace and plot streamlines
numStreamlines = 500;
clf
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)-scrsz(3)/3 scrsz(3) scrsz(3)/3]);
subplot(1,3,1);
title('Streamlines');
h = plotGrid(G, 'facea', 0.3, 'edgea',0.1);
axis equal tight;
pos = [ones(numStreamlines,1), rand(numStreamlines,2)];
timer=tic;
[xyz,t,c] = pollock(G, x, pos, 'substeps', 1, 'maxsteps',10000);
thetime=toc(timer);
dispif(true, 'Traced %d steps in %d streamlines in %g second\n', ...
   sum(cellfun(@numel, t)), numStreamlines, thetime);
streamline(xyz);

%% Extract time-of-flight
subplot(1,3,2);
title('Streamline time of flight')

% add time of flight along streamlines
ct = cellfun(@cumsum, t, 'uniformoutput', false);

% rearrange to plain double arrays
t  = vertcat(t{:});
i= (t ~= inf);
t  = t(i);
ct = vertcat(ct{:}); ct = ct(i);
c  = vertcat(c{:});  c  = c(i);

% Compute weigthed average of time of flight in each cell that has been
% crossed by one or more streamlines,
%
%    tof_i = sum(t_ij*ct_ij; j) / sum(t_ij; j)
%
% where t_ij is the time-of-flight diff through cell i of streamline j

tof = accumarray(c, ct.*t, [G.cells.num, 1], [], inf)./...
   accumarray(c, t,     [G.cells.num, 1], [], 1);


plotCellData(G, tof, 'EdgeColor','none');
axis equal tight;

%% Compute time-of-flight using diagnostics tools
subplot(1,3,3);
title('Finite-volume time of flight')
fdTOF = computeTimeOfFlight(x, G, rock, 'src', src);
plotCellData(G, fdTOF,'EdgeColor','none');
axis equal tight;

% Use fdTOF for caxis scalling as this does not change randomly.
subplot(1,3,2), caxis([0,fdTOF(src.cell(end))])
subplot(1,3,3), caxis([0,fdTOF(src.cell(end))])

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
