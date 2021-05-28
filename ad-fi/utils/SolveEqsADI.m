function [dx, linsolver_diverged] = SolveEqsADI(eqs, phi)
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

useBasis = ~isempty(phi);

if useBasis
    % We have been provided a pod basis, use it to reduce equation sets
    % before solving.
    for i = 1:numel(eqs)
        eqs{i}.val = phi.basis{i}'*eqs{i}.val;
    end
end

numVars = cellfun(@numelValue, eqs)';
cumVars = cumsum(numVars);
ii = [[1;cumVars(1:end-1)+1], cumVars];

eqs = cat(eqs{:});
tic

% Above CAT means '.jac' is a single element cell array.  Extract contents.
J = -eqs.jac{1};

if useBasis
    blkphi = blkdiag(phi.basis{:});
    J = blkphi'*J*blkphi;
end

tmp = J\eqs.val;

linsolver_diverged = false;
if ~all(isfinite(tmp))
   linsolver_diverged = true;
   warning('Linear solver produced non-finite values! Stop simulation.\n');
   dx = [];
   return
end
   
eqn = size(ii,1);
dx = cell(eqn,1);
for i = 1:eqn
    dx{i} = tmp(ii(i,1):ii(i,2));
end

if useBasis
    for i = 1:numel(dx)
        dx{i} = phi.basis{i}*dx{i};
    end
end
