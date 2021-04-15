function [z, dz] = testFunction(x,y, varargin)
% Simple quadratic function:
%      z(x,y) = ((x-x0)/rx)^2 + ((y-y0)/ry)^2
% possibly rotated by angle th around (x0,y0)
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
opt = struct('x0', .5, ...
             'y0', .5, ...
             'rx',  1, ...
             'ry',  1, ...
             'th',  0);
opt = merge_options(opt, varargin{:});

outputDerivatives = nargout > 1;
if outputDerivatives
    [x,y] = initVariablesADI(x,y);
end

[x0, y0] = deal(opt.x0, opt.y0);
[rx, ry] = deal(opt.rx, opt.ry);
th = opt.th;



xx =  cos(th)*(x-x0) + sin(th)*(y-y0);
yy = -sin(th)*(x-x0) + cos(th)*(y-y0);

z  = -(xx/rx).^2 - (yy/ry).^2;
if outputDerivatives
    dz = [z.jac{1}, z.jac{2}]';
    z  = z.val;
end
%z  = -((x-x0)/rx).^2 - ((y-y0)/ry).^2 + .25*cos((x-x0)/rx+y-y0);
%dz = -[2*(x-x0)/rx^2+.25*sin((x-x0)/rx+y-y0), 2*(y-y0)/ry^2+.25*sin((x-x0)/rx+y-y0)]';
end

