function [qSurf, BCTocellMap, cells] = getSourceFluxesAD(model, mob, s, src)
%Short description
%
% SYNOPSIS:
%   [qSurf, BCTocellMap, cells] = getSourceFluxesAD(model, mob, b, s, src)
%
% DESCRIPTION:
%   
%
% REQUIRED PARAMETERS:
%   model          - Subclass of ReservoirModel indicating which phases are
%                    active.
%
%   mob            - A cell array of cell mobility values for all active
%                    phases.
%
%
%   s              - A cell array of saturations per cell for all active
%                    phases.
%
%   src            - Source struct as defined by addSource
%
% RETURNS:
%   qSurf           - Source terms at standard conditions. Cell array of
%                   same dimensions as the number of active phases.
%
%   cells          - A list of cells for which the entries of qRes should
%                    be added.

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

    nPh = sum(model.getActivePhases);
    assert(size(src.sat, 2) == nPh);
    
    
    cells = src.cell;
    nsrc = numel(cells);
    
    inj = src.rate > 0;
    prod = ~inj;
    qSurf = cell(nPh, 1);
    
    if any(prod)
        totMob = 0;
        for i = 1:nPh
            totMob = totMob + mob{i};
        end
    end
    
    for i = 1:nPh
        q = zeros(nsrc, 1);
        if any(prod)
            q = double2ADI(q, totMob);
        end
        
        if any(inj)
            % Injection rates are given in reservoir conditions
            q(inj) = src.rate(inj).*src.sat(inj, i);
        end
        
        if any(~inj)
            c = cells(~inj);
            % Use mobilities ratios to ensure that immobile fluids are not
            % removed from the reservoir.
            f = mob{i}(c)./totMob(c);
            q(~inj) = q(~inj) + f.*src.rate(~inj);
        end
        qSurf{i} = q;
    end
    cellToBCMap = sparse((1:nsrc)', cells, 1, nsrc, model.G.cells.num);
    BCTocellMap = cellToBCMap';
end
