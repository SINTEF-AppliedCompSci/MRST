%% Small test for the simple equilibrium routine
G = cartGrid([10, 1, 10], [1, 1, 1]);
G = computeGeometry(G);

% First contact at .37, second at .8 for a total of three phases present
contacts = [.37, .8];

% Compute equilibrium saturations
s = simpleEquilibrium(G, contacts);

% Plot the different saturations, as well as the different contacts
figure(1); clf
for i = 1:size(s, 2)
    subplot(1, size(s, 2), i);
    plotCellData(G, s(:, i))
    caxis([0, 1])
    view(0, 0)
    hold on
    for j = 1:numel(contacts)
        plot3([0, 1], [-0.01, -0.01], contacts([j, j]), 'r', 'linewidth', 3)
    end
end

%% The equilibriation does not have to happen along the z-axis
% By specifying a vector, we can define the saturations to be aligned with
% the force of gravity in any direction.
vec = [1, 0, 1];
s = simpleEquilibrium(G, contacts, vec);

figure(1); clf
for i = 1:size(s, 2)
    subplot(1, size(s, 2), i);
    plotCellData(G, s(:, i))
    caxis([0, 1])
    view(0, 0)
    hold on
end

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
