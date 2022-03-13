function vi = velocityInterpolation(G, type)
    % Construct velocity vector from face fluxes

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    switch type
        case 'mimetic'

            cellNo    = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
            faceNo = G.cells.faces(:,1);
            C = G.faces.centroids(faceNo, :) - G.cells.centroids(cellNo, :);
            sgn = 1 - 2*(cellNo ~= G.faces.neighbors(faceNo,1));

            D = cell(1,G.griddim);
            for dNo = 1:G.griddim
                D{dNo} = sparse(cellNo, faceNo, C(:,dNo).*sgn, G.cells.num, G.faces.num);
            end

            f2c = @(v) faceFlux2cellVelocity(D,v);
            vi = struct('D', {D}, 'faceFlux2cellVelocity', f2c);
            
        otherwise
                error('Unknown velocity interpolation type')
    end

end

function vc = faceFlux2cellVelocity(D, v)
    
    vc = cell(1,numel(D));
    for dNo = 1:numel(D)
        vc{dNo} = D{dNo}*v;
    end
    vc = SpatialVector(vc{:});

end
