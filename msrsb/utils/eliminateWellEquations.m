function [A_pp, q_p, A_ww, A_wp, q_w] = eliminateWellEquations(A, q, nc)
% Eliminate well equations from linear system

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
    iscell = false(size(q));
    iscell(1:nc) = true;
    
    A_pp = A(iscell, iscell);
    A_pw = A(iscell, ~iscell);
    q_p = q(iscell);
    
    A_ww = A(~iscell, ~iscell);
    A_wp = A(~iscell, iscell);
    q_w = q(~iscell);
    
    A_pp = A_pp - A_pw*(A_ww\A_wp);
    q_p = q_p - A_pw*(A_ww\q_w);
end
