function T = getFaceTransmissibility(G, rock, deck, varargin)
%Compute face transmissibilities, accounting for input-specific multipliers
%
% SYNOPSIS:
%    T = getFaceTransmissibility(G, rock, deck)
%    T = getFaceTransmissibility(G, rock)
%
% DESCRIPTION:
%   Computes transmissibilities per interface. Accounts for multipliers due
%   to both general transmissibility multipliers and fault multipliers
%   specifically.
%
% REQUIRED PARAMETERS:
%   G    - Valid grid structure.
%
%   rock - Valid rock structure.
%
%   deck - (OPTIONAL) ECLIPSE style input deck used to produce the grid,
%                     typically produced by readEclipseDeck. Needed to
%                     account for multipliers.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   (Passed directly onto underlying function computeTrans)
%
% RETURNS:
%   T    - Transmissibilities, one per interface.
%
% SEE ALSO:
%   computeTrans

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
    if nargin == 2 || isempty(deck)
        deck = struct('GRID', struct());
    end
    
    if mod(numel(varargin), 2) == 1
        varargin = {deck, varargin{:}};
        deck = struct('GRID', struct());
    end
    T = computeTrans(G, rock, varargin{:});
    
    % Multiply inn transmissibility multipliers
    m = computeTranMult(G, deck.GRID);
    if ~isempty(m)
        T = m.*T;
    end
    % Reduce half-face trans for face transmissibility
    cf = G.cells.faces(:,1);
    nf = G.faces.num;
    T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
    
    % Treat fault multipliers
    flt = processFaults(G, deck.GRID);
    if ~isempty(flt)
        for i = 1:numel(flt)
            ff = flt(i).faces;
            % Max of zero is to avoid negative values or NaNs
            T(ff) = T(ff).*max(flt(i).mult, 0);
        end
    end
end
