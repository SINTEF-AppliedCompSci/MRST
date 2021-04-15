function [A, qPluss, signQ] = ...
        generateUpstreamTransportMatrix(G, S, W, resSol, wellSol, varargin)
% generateUpstreamTransportMatrix for use in saturation solver
%
% SYNOPSIS:
% [A, qPluss, signQ] = generateUpstreamTransportMatrix(G, S, W, resSol, ...
%                                                    wellSol, pn1, pv1, ...)
%
% DESCRIPTION:
%   Generates sparse matrix A(cellFlux) which is used in transport solver
%   s1 = s0 + dt*Dv*(Af(s) + q+).
%   Assumes no-flow boundary conditions on all cell-faces
%
%   signQ is useful for differentation of max(q, 0)/min(q, 0) wrt to q,
%   when q is zero
%
% REQUIRED PARAMETERS:
%
% OPTIONAL PARAMETERS:
%   Transpose     - if true, the transpose of a is given (default false)
%
%   VectorOutput  - if true, output A is a struct with fields 'i', 'j' and
%                   'qMinus', such that
%                   A = sparse(i, j, -cellFlux) + diag(qMinus)
%                   (NOTE MINUSES). Default value is false
%
%   RelativeThreshold  - Considers values below
%                        max(abs(cellflux))*RealtiveThreshold
%                        as zero (default 0)

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


% ---------------------------------------------
opt = struct('Transpose'         , false, ...
             'VectorOutput'      , false, ...
             'RelativeThreshold' , 0);
opt = merge_options(opt, varargin{:});

% ---------------------------------------------

cellFlux    = faceFlux2cellFlux(G, resSol.flux);
wellRates   = vertcat(wellSol.flux);
wellSys     = [ W.S ];
Cw          = vertcat(wellSys.C);

cfns  = cellFaceNeighbors(G);
ii    = cfns(:,1);
q           = Cw'*wellRates;

if opt.RelativeThreshold > 0
    maxV              = max( abs(cellFlux) );
    belowThresholdInx = ( abs(cellFlux) < opt.RelativeThreshold/maxV );
    cellFlux(belowThresholdInx)  = 0;
    belowThresholdInx = ( abs(q) < opt.RelativeThreshold/maxV );
    q(belowThresholdInx)         = 0;
end

negFlux     = and( cellFlux < 0, cfns(:,2));  % only allow interior faces

%extInd      = (ii == 0);
%ii(extInd)  = cfns(extInd, 2);

jj          = ii;
jj(negFlux) = cfns(negFlux, 2);

qMinus        = q;
qMinus(q > 0) = 0;

% -----------------------------------

if opt.VectorOutput
    A.i      = ii;
    A.j      = jj;
    A.qMinus = qMinus;
else
    numC  = G.cells.num;
    if opt.Transpose
        A = sparse(jj, ii, -cellFlux, numC, numC) + ...
            spdiags(qMinus, 0, numC, numC);
    else
        A = sparse(ii, jj, -cellFlux, numC, numC) + ...
            spdiags(qMinus, 0, numC, numC);
    end
end

if nargout > 1
    qPluss = q - qMinus;
end

%{
if nargout > 2
    Dw    = blkdiag( wellSys.D );
    if isfield(W, 'sign')
        signs = vertcat(W.sign);
    else
        signs    = ones( numel(W), 1 );
        totRates = Dw'*wellRates;
        signs( totRates < 0 ) = -1;
    end
    signQ       = Cw'*Dw*signs;
end
%}
if nargout > 2
    Dw    = blkdiag( wellSys.D );
    signs    = ones( numel(W), 1 );
    totRates = Dw'*wellRates;
    signs( totRates < 0 ) = -1;
    for k = 1 : numel(W)
        if ~isempty(W(k).sign) % override if sign is given expl
            signs(k) = W(k).sign;
        end
    end
    signQ       = Cw'*Dw*signs;
end

function cfns = cellFaceNeighbors(G)
% finds neighbor pairs according to cellface - cells
% [cellface-cells celface-cell-neigbors
numFaces = diff(G.cells.facePos);

cfns        = double( G.faces.neighbors(G.cells.faces(:,1), :) );
cellNo      = rldecode(1:G.cells.num, double(numFaces), 2)';
sgn         = 2*(G.faces.neighbors(G.cells.faces(:,1), 1) == cellNo)-1;
flipRowInx  = ( sgn < 0 );
cfns(flipRowInx, :) = cfns(flipRowInx, [2 1]);
extInd              = find(cfns(:,1)==0);
cfns(extInd, :)     = cfns(extInd, [2 1]);
return

