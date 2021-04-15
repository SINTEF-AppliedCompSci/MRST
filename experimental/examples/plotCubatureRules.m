close all;
k = 13;

%%

[x, w, n, xR] = getLineCubaturePointsAndWeights(k);

figure('Name', 'Line cubature')
hold on
line([-1,1], [0,0]);
wMax = max(w);
wMin = min(w);
nClr = 40;
clr = jet(nClr);
cw = @(w) 50/(wMax - wMin)*(w - wMin) + 10;
for pNo = 1:numel(w)
    cIx = floor((w(pNo) - wMin)./(wMax - wMin)*(nClr-1)) + 1;
    plot(x(pNo,1), 0, '.', 'markerSize', cw(w(pNo)), 'color', clr(cIx,:));
end
axis equal tight

%%

[x, w, n, xR] = getTriangleCubaturePointsAndWeights(k);

figure('Name', 'Triangle cubature')
hold on
patch('faces', [1,2,3], 'vertices', xR, 'facec', 'none');
wMax = max(w);
wMin = min(w);
nClr = 40;
clr = jet(nClr);
cw = @(w) 50/(wMax - wMin)*(w - wMin) + 10;
for pNo = 1:numel(w)
    cIx = floor((w(pNo) - wMin)./(wMax - wMin)*(nClr-1)) + 1;
    plot(x(pNo,1), x(pNo,2), '.', 'markerSize', cw(w(pNo)), 'color', clr(cIx,:));
end
axis equal tight

%%

[x, w, n, xR] = getTetrahedronCubaturePointsAndWeights(k);

figure('Name', 'Tetrahedron cubature')
hold on
patch('faces', [1,2,3; 1,2,4; 1,3,4; 2,3,4], 'vertices', xR, 'facec', 'none');

wMax = max(w);
wMin = min(w);
nClr = 40;
clr = jet(nClr);
cw = @(w) 50/(wMax - wMin)*(w - wMin) + 10;
for pNo = 1:numel(w)
    cIx = floor((w(pNo) - wMin)./(wMax - wMin)*(nClr-1)) + 1;
    plot3(x(pNo,1), x(pNo,2), x(pNo,3), '.', 'markerSize', cw(w(pNo)), 'color', clr(cIx,:));
end
view([120,30]); axis equal tight

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
