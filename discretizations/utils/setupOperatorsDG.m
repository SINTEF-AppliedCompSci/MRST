function op = setupOperatorsDG(disc, G, rock, varargin)
%Undocumented Utility Function

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

    op = setupOperatorsTPFA(G, rock, varargin{:});
    
    [~, ~, ~, f] = disc.getCubature(find(disc.internalConn), 'face'); %#ok
    nf = numel(f);
    
    N = [1:nf; nf+1:2*nf]';
    
    M = sparse((1:nf)'.*[1,1], N, 0.5, nf, 2*nf);
    op.M = M;
    op.faceAvg = @(x) M*x;
    
    op.faceUpstr0 = @(flag, x) op.faceUpstr(flag, x);
    op.faceUpstr = @(flag, x) faceUpstr(flag, x, N, [nf, 2*nf]);
    op.faceUpstr = @(flag, x) fup(op, flag, x, nf);
    
    op.C0 = op.C;
    C = sparse((1:nf)'.*[1,1], N, ones(nf,1)*[1 -1], nf, 2*nf);
    op.C = C;
    
    op.Grad = @(x) grad(op, x);
    
    op.Div = @(v) div(op, v);
    
    op.T = op.T_all(f);
    
    op.velocityInterp = velocityInterpolation(G, 'mimetic');
    
end

function dx = grad(op, x)
    if size(value(x),1) < size(op.C,2)
        dx = -op.C0*x;
    else
        dx = -op.C*x;
    end
end

function dv = div(op, v)
    if size(value(v),1) < size(op.C',2)
        dv = op.C0'*x;
    else
        dv = op.C'*x;
    end
end

function up = fup(op, flag, x, nf)
    if numel(flag) < nf
        up = op.faceUpstr0(flag, x);
    else
        up = op.faceUpstr(flag, x);
    end
end
