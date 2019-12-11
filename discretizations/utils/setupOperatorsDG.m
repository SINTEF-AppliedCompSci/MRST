function op = setupOperatorsDG(disc, G, rock, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

    op = setupOperatorsTPFA(G, rock, varargin{:});
    
    [~, ~, ~, f] = disc.getCubature(find(disc.internalConn), 'face');
    nf = numel(f);
    
    N = [1:nf; nf+1:2*nf]';
    
    M = sparse((1:nf)'.*[1,1], N, 0.5, nf, 2*nf);
    op.M = M;
    op.faceAvg = @(x) M*x;
    
    op.faceUpstr = @(flag, x) faceUpstr(flag, x, N, [nf, 2*nf]);
    
    C = sparse((1:nf)'.*[1,1], N, ones(nf,1)*[1 -1], nf, 2*nf);
    op.C = C;
    op.Grad = @(x) -C*x;
    
    op.Div = @(x) -C*x;
    
    op.T = op.T_all(f);
    
    op.velocityInterp = velocityInterpolation(G, 'mimetic');
    
end
