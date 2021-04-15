function p = msfvPressureSolve(arg1, rhs, DG, useCorrection)
%Undocumented Utility Function

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

    persistent B Cx updateTargets clusters A_nn;
    if nargin == 0
        % Reset all basis functions
        B = [];
        A_nn = [];
        Cx = [];
        return
    end

    if nargin == 1
        % Flag cells for update before next solve
        updateTargets = arg1;
        return
    end

    A = arg1;
    nc = size(A, 1);

    if isempty(B)
        if useCorrection && isempty(Cx)
            [B, Cx] = createMSFVBasis(A, DG, useCorrection);
        else
            B = createMSFVBasis(A, DG, useCorrection);
        end
        updateTargets = false(nc, 1);
    end

%     if any(updateTargets)
%         % Do update
%         if isempty(clusters)
%
%             clusters = zeros(nc, 1);
%             ni = numel(DG.ii);
%
%             Ap = DG.P*A*DG.P';
%             Ap = Ap(1:ni, 1:ni);
%             clusters(DG.ii) = components(Ap);
%         end
%         B = createMSFVBasis(A, DG, updateTargets, clusters);
%         updateTargets = false(size(A, 1), 1);
%     end
    X = DG.X;

    if isempty(A_nn)
        A_nn = X*A*B;
    end
    if useCorrection
        C_rhs = Cx(rhs);
        q = X*(rhs - A*C_rhs);
    else
        C_rhs = zeros(size(rhs));
    end

    p = B*mldivide(A_nn, q) + C_rhs;

end
