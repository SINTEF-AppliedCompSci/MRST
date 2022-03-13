classdef MomentFitting3DCubature < Cubature
    % Cubature based on moment-fitting for MRST grids
    
    properties
        
        reduce
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = MomentFitting3DCubature(G, prescision)
            % Set up cubatureature
            
            % Basic properties handled by parent class
            cubature = cubature@Cubature(G, prescision, 3);
        end
           
        %-----------------------------------------------------------------%
        function cubature = makeCubature(cubature)
            
            % Dimension of cubature
            dim = 3;
            G   = cubature.G;
            % Basis functions used in moment-fitting
            basis = dgBasis(dim, cubature.prescision, 'legendre');
            psi   = basis.psi;
            nDof  = basis.nDof;

            if cubature.prescision <=1
                % Precision is 1, use midpoint rule
                [x, ~, n] = getCubeCubaturePointsAndWeights(cubature.prescision);
                w = G.cells.volumes;
                x = repmat(x, G.cells.num,1);
                n = repmat(n, G.cells.num,1);
                cellNo = (1:G.cells.num)';
                
            else
                % If precision is higher than one, we must calculate
                % weights based on moment-fitting
                
                % The starting point is a quadrature for a reference square
                % with more than nDof points
                if 1
                    G1 = computeGeometry(cartGrid([1,1,1], [2,2,2]));
                    G1.nodes.coords = G1.nodes.coords - 1;
                    G1 = computeVEMGeometry(G1);
                    G1 = computeCellDimensions(G1);
                    cubTet = TetrahedronCubature(G1, cubature.prescision);
                    [~, x, ~, cellNo, ~] = cubTet.getCubature(1, 'cell');
                    x = cubTet.transformCoords(x, cellNo);
                    x = unique(x, 'rows');
                 else
                    k = 0;
                    nS = 0;
                    while nS < 2*nDof
                        [x, ~, nS] = getCubeCubaturePointsAndWeights(cubature.prescision + k);
                        k = k+1;
                    end
                end
                
                % We use known cubature to calculate the moments
                if isfield(G, 'parent')
                    knownCub = CoarseGrid3DCubature(G, cubature.prescision);
                else
                    knownCub = TetrahedronCubature(G, cubature.prescision);
                end
                [~, xq, wq, cellNo] = knownCub.getCubature((1:G.cells.num)', 'cell');
                xq     = cubature.transformCoords(xq, cellNo);
                % Moments
                M = cellfun(@(p) accumarray(cellNo, wq.*p(xq)), psi, 'unif', false);
                % Compute right-hand side
                rhs = zeros(nDof, G.cells.num);
                tol = eps(1);
                for dofNo = 1:nDof
                    m = M{dofNo};
                    m(abs(m) < tol) = 0;
                    rhs(dofNo, :) = m;
                    M{dofNo} = m;
                end
                [x,w,n] = fitMoments(x, basis, M, 'equal', G.cells.equal);
                cellNo = rldecode((1:G.cells.num)', n, 1);
                w = w.*G.cells.volumes(cellNo);
            end
            
            % Map from reference to physical coordinates
            x = cubature.transformCoords(x, cellNo, true);
            cubature = cubature.assignCubaturePointsAndWeights(x,w,n);
            
        end        
    end
    
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
