function fn = getUniformInterpRegLinear(xs, fs, reg)
% Intermediate interpolator. Support for multiple regions, with the caveat
% that all functions must be given on the same set of points.

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

    ns = numel(xs);
    assert(issorted(xs));
    assert(size(xs, 1) == 1);
    assert(size(xs, 2) == size(fs, 2));
    assert(max(reg) <= size(fs, 1));
    assert(all(reg) > 0);
    
    slope = bsxfun(@rdivide, diff(fs, 1, 2), diff(xs, 1, 2));
    slope = slope(:, [1:(ns-1), ns-1]);
    
    % Dynamic part
    fs_unpack = reshape(fs', 1, [])';
    slope_unpack = reshape(slope', 1, [])';
    fn = @(x) interpReg(x, xs, fs_unpack, slope_unpack, reg);
end

function f = interpReg(x, xs, fs, slope, reg)
    xv = double(x);

    nx = numel(xv);
    ns = numel(xs);
    
    act = bsxfun(@le, xv, xs);
    [sel, jj] = max(act, [], 2);
    jj(~sel) = ns;
    jj = max(jj - 1, 1);

    
    nreg = size(fs, 1)/ns;
    % Unpack regions, since the grid is assumed uniform
    pick = sparse((1:nx)' , jj + (reg-1)*ns, 1, nx, ns*nreg);
    f0 = pick*fs;
    dfdx = pick*slope;

    f =  f0 + dfdx.*(x - xs(jj)');
end
