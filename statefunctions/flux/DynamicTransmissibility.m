classdef DynamicTransmissibility < StateFunction
%State function for computing dynamic transmissibilities

    properties
        harmonicAvgOperator
        twoPointOperator
        conductivity_name
    end
    
    methods
        %-----------------------------------------------------------------%
        function dt = DynamicTransmissibility(model, conductivity_name)
            if nargin < 2
                conductivity_name = 'Permeability';
            end
            dt@StateFunction(model);
            dt = dt.dependsOn(conductivity_name);
            dt.conductivity_name   = conductivity_name;
            dt.twoPointOperator    = getTwoPointOperator(model.G);
            dt.harmonicAvgOperator = getHarmonicAvgOpeartor(model.G);
            dt.label = 'T';
        end
        
        %-----------------------------------------------------------------%
        function T = evaluateOnDomain(prop, model, state, allFaces)
            lambda = prop.getEvaluatedDependencies(state, prop.conductivity_name);
            if iscell(lambda)
                T = cell(numel(lambda),1);
                for i = 1:numel(lambda)
                    if ~isempty(lambda{i})
                        T{i} = prop.getTransmissibility(lambda{i});
                    end
                end
            else
                T = prop.getTransmissibility(lambda);
                if (nargin < 4 || ~allFaces)
                    T = T(model.operators.internalConn);
                end
            end
        end
        
        %-----------------------------------------------------------------%
        function T = getTransmissibility(prop, lambda)
            % Compute two-point difference and to harmonic average
            T = prop.twoPointOperator(lambda);
            T = prop.harmonicAvgOperator(T);
            % Handle negative transmissibility
            fix    = T < 0;
            T(fix) = -T(fix);
        end
    end
end

%-------------------------------------------------------------------------%
function tp = getTwoPointOperator(G)
    % Mappings from cells to its faces
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces = G.cells.faces(:,1);
    % Vector from cell to face centroid
    C = G.faces.centroids(faces,:) - G.cells.centroids(cells,:);
    % Oriented normals
    sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
    N   = bsxfun(@times, sgn, G.faces.normals(faces, :));
    % Make function
    cn  = sum(C.*N,2)./sum(C.*C,2);
    tp = @(lambda) cn.*lambda(cells);
end

%-------------------------------------------------------------------------%
function ha = getHarmonicAvgOpeartor(G)
    % Harmonig averaging operator
    faces = G.cells.faces(:,1);
    M = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
    ha = @(T) 1./(M*(1./T));
end

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