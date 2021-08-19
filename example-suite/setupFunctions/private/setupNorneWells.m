function [W] = setupNorneWells(G,rock)
% Helper function to generate a simplified well setup for the Norne model.
%
% SYNOPSIS:
%   [W] = setupNorneWells(G,rock))
%
% DESCRIPTION:
%   
%    The reservoir is produced using a set of 6 production wells controlled by
%    bottom-hole pressure and 6 rate-controlled injectors. 
%    For simplicity, all wells are assumed to be vertical and are
%    assigned using the logical (i,j) subindex.
%    Note this well setup is different to the well setup provided in
%    NORNE.GRDECL (as used in showNorne.m).
%
% PARAMETERS:
%    G - grid structure for Norne generated from setupNorneRealization.
%    rock - rock structure for Norne generated from setupNorneRealization.
%
% RETURNS:
%    W - well structure for G and rock.
%
% SEE ALSO:
%   `setupNorneRealization`, `setupNorneFn`
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


% Set vertical injectors, completed in the lowest 12 layers.
nz = G.cartDims(3);
I = [ 9, 26,  8, 25, 35, 10];
J = [14, 14, 35, 35, 68, 75];
R = 4*[ 4,  4,  4,  4,  4,  4]*1000*meter^3/day;
P = [300, 300, 300, 300, 300, 300]*barsa;
W = [];
for i = 1 : numel(I)
    W = verticalWell(W, G, rock, I(i), J(i), nz-11:nz, 'Type', 'rate', ...
        'InnerProduct', 'ip_tpf', 'Sign', 1, ...
        'Val', R(i), 'Radius', 0.1, 'Comp_i', [1, 0], ...
        'name', ['I', int2str(i)], 'refDepth', 2500);
end

% Set vertical producers, completed in the upper 14 layers
I = [17, 12, 25, 35, 15];
J = [23, 51, 51, 95, 94];

R = -sum(R)/numel(I)*ones(numel(I),1);
P = 100*barsa*ones(numel(I),1);

for i = 1 : numel(I)
    W = verticalWell(W, G, rock, I(i), J(i), 1:14, 'Type', 'bhp', ...
        'InnerProduct', 'ip_tpf', 'Sign', -1, ...
        'Val', P(i), 'Radius', 0.1, 'refDepth', 2500, ...
        'name', ['P', int2str(i)], 'Comp_i', [0, 1]);
end

