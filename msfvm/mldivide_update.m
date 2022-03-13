function x = mldivide_update(A, b, x_old, targets, complist)
% A,b linear system
% A_old, b_old old linear system
% x_old old solution vector
% targets, list of cells which should be updated (linear subsystems
% containing one of these cells must be updated)
% complist as defined by components(A)
% A should have the same graph, but not the same values as A_old

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


target_blocks = complist(targets);
ind = ismember(complist, target_blocks);
x_old(ind,:) = mldivide(A(ind,ind), b(ind,:));
x = x_old;
end
