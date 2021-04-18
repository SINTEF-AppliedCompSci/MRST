% Test marking cut domains from sliceGrid

clear
close all

legacy = verLessThan('matlab', '8.6');

nx = 5;
physdims = [1, 1, 1];
G0 = cartGrid([nx, nx, nx], physdims);
G0 = computeGeometry(G0);
n = [1, 1, 0];

%% Single cut
[G, gix] = sliceGrid(G0, physdims/2, 'normal', n);
start = 0.1*physdims;
if legacy
    m = markCutGrids(G, gix.new.faces, 'legacy', true, 'start', [0.1,0.1,0.1]);
else
    m = markCutGrids(G, gix.new.faces);
end
figure
plotCellData(G, m);
view(3)

%% Two parallel cuts, three domains
offset = 0.1;
cuts = [physdims/2; physdims/2+offset];
[G, gix] = sliceGrid(G0, cuts, 'normal', n);
if legacy
    x0 = [0.1, 0.1, 0.1];
    m1 = markCutGrids(G, gix.new.faces, 'legacy', true, 'start', 1-x0);
    m2 = markCutGrids(G, gix.new.faces, 'legacy', true, 'start', x0);
    m = m1+2*m2;
else
    m = markCutGrids(G, gix.new.faces);
end
figure
plotCellData(G, m);
view(3)

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
