mrstModule add vem mpfa mpsaw vemmech libgeometry

close all

%% params setting
% nref     : degree of refinement
% Nd       : dimension (2D or 3D)
% kappa    : Value of heterogenity in the domain (see paper)
% alpha    : Parameter to setup Lamé coefficient lambda (lambda = alpha*mu)
% gridtype : grid type (see mpsaPaperConvergenceFunc)
% eta      : Value used to set the position of the continuity point

% Possibility to run vem for comparison
doVem = false;

%% New Case
dothiscase = true;
params = struct('nref'    , 4, ...
                'Nd'      , 3, ...
                'kappa'   , 10, ...
                'alpha'   , 0, ...
                'gridtype', 5, ... % Cartesian
                'eta'     , 1/4);

output = mpsaPaperConvergenceFunc(params, 'doVem', doVem, 'blocksize', 100);

figure
hold on
plotConvTest(output, params);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% The MPSA-W module is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% The MPSA-W module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with the MPSA-W module.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>

