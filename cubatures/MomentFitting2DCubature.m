classdef MomentFitting2DCubature < Cubature
    % Cubature based on moment-fitting for MRST grids
    
    properties
        reduce = true;
        chunkSize = 10;
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function cubature = MomentFitting2DCubature(G, prescision, varargin)
            % Set up cubature
            
            % Most of the construction is handled by parent class
            cubature        = cubature@Cubature(G, prescision, 2, varargin{:});
            cubature        = merge_options(cubature, varargin{:});
        end
           
        %-----------------------------------------------------------------%
        function cubature = makeCubature(cubature)
            
            % Dimension of cubature
            dim = 2;
            G   = cubature.G;
            % Basis functions used in moment-fitting
            basis = dgBasis(dim, cubature.prescision, 'legendre');
            psi   = basis.psi;
            nDof  = basis.nDof;
            
            if cubature.prescision <= 1
                % Precision is 1, use midpoint rule
                [x, ~, n] = getSquareCubaturePointsAndWeights(cubature.prescision);
                if G.griddim == 2
                    w = ones(G.cells.num,1);
                    x = repmat(x, G.cells.num, 1);
                    n = repmat(n, G.cells.num, 1);
                    type = 'cell';
                else
                    w = ones(G.faces.num,1);
                    x = repmat(x, G.faces.num, 1);
                    n = repmat(n, G.faces.num, 1);
                    type = 'face';
                end
            else
                % If precision is higher than one, we must calculate
                % weights based on moment-fitting
                
                % The starting point is a quadrature for a reference square
                % with more than nDof points
                if 1
                    G1 = computeGeometry(cartGrid([1,1], [2,2]));
                    G1.nodes.coords = G1.nodes.coords - 1;
                    G1 = computeVEMGeometry(G1);
                    G1 = computeCellDimensions(G1);
                    cubTri = TriangleCubature(G1, cubature.prescision+1);
                    [~, x, ~, cellNo, ~] = cubTri.getCubature(1, 'cell');
                    x = cubTri.transformCoords(x, cellNo);
                    x = unique(x, 'rows');
                else
                    k = 0;
                    nS = 0;
                    while nS < 2*nDof
                        [x, ~, nS] = getSquareCubaturePointsAndWeights(cubature.prescision + k);
                        k = k+1;
                    end
                end
                % Cubature is either for faces or cells
                if G.griddim > cubature.dim
                    type = 'face';
                    elements = 1:G.faces.num;
                    equal = false;
                    if isfield(G.faces, 'equal')
                        equal = G.faces.equal;
                    end
                else
                    type = 'cell';
                    elements = 1:G.cells.num;
                    equal = false;
                    if isfield(G.cells, 'equal')
                        equal = G.cells.equal;
                    end
                end
                % We use known cubature to calculate the moments
                if isfield(G, 'parent')
                    knownCub = CoarseGrid2DCubature(G, cubature.prescision);
                else
                    knownCub = TriangleCubature(G, cubature.prescision);
                end
                [~, xq, wq, cellNo, faceNo] = knownCub.getCubature(elements, type);
                % Map cubature points to reference coordinates
                if G.griddim == 3
                    % Map to face reference coordinates
                    xq = G.faces.phys2ref(xq, faceNo);
                    count = faceNo;
                else
                    % Map to cell reference coordiantes
                    xq    = cubature.transformCoords(xq, cellNo);
                    count = cellNo;
                end
                % Moments
                tol = eps(10);
                M = cellfun(@(p) accumarray(count, wq.*p(xq)), psi, 'unif', false);
                M = cellfun(@(m) m.*(abs(m) > tol), M, 'unif', false);
                [x,w,n] = fitMoments(x, basis, M, 'equal', equal, 'reduce', cubature.reduce, 'chunkSize', cubature.chunkSize);
            end
            
            % Map from reference to physical coordinates
            if strcmp(type, 'face')
                % Face coordinates
                faceNo = rldecode((1:G.faces.num)', n, 1);
                x = G.faces.ref2phys(x, faceNo);
                w = w.*G.faces.areas(faceNo);
            else
                % Cell coordinates
                cellNo = rldecode((1:G.cells.num)', n, 1);
                x = cubature.transformCoords(x, cellNo, true);
                w = w.*G.cells.volumes(cellNo);
            end
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
