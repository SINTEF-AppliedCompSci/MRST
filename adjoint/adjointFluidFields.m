function fluid = adjointFluidFields(fluid)
%Extend fluid functionality with fields needed in (2nd order) adjoint imp.
%
% SYNOPSIS:
%   fluid = adjointFluidFields(fluid)
%
% PARAMETERS:
%   fluid - A 'fluid_structure'.
%
% RETURNS:
%   fluid - Updated fluid structure now containing additional fields
%             - dLtInv (state) -- d/ds     (1 / \lambda_t(state))
%             - d2LtInv(state) -- d^2/ds^2 (1 / \lambda_t(state))
%             - d2fw   (state) -- d^2/ds^2 f_w(state)
%
% SEE ALSO:
%   `fluid_structure`.

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


   fluid.dLtInv  = @(state) dLtInv (state, fluid);
   fluid.d2LtInv = @(state) d2LtInv(state, fluid);
   fluid.d2fw    = @(state) d2fw   (state, fluid);
end

%--------------------------------------------------------------------------

function r = dLtInv(state, fluid)
% (d/ds) (1 / \lambda_t(s))
   [mob, dmob] = mobilities(state, fluid);

   Lt =   sum(mob,  2);
   r  = - sum(dmob, 2) ./ Lt.^2;
end

%--------------------------------------------------------------------------

function r = d2LtInv(state, fluid)
% (d^2/ds^2) (1 / \lambda_t(s))
   [mob{1:3}] = mobilities(state, fluid);

   Lt  = sum(mob{1}, 2);
   dLt = sum(mob{2}, 2);

   r = (2./Lt .* dLt.^2 - sum(mob{3}, 2)) ./ Lt.^2;
end

%--------------------------------------------------------------------------

function d2f = d2fw(state, fluid)
% (d^2/ds^2) f(s)
   [mob{1:3}] = mobilities(state, fluid);

   Lt  = sum(mob{1}, 2);
   f   = bsxfun(@rdivide, mob, Lt);
   df  = (f(:,2).*mob{2}(:,1) - f(:,1).*mob{2}(:,2)) ./ Lt;   % Chain rule.
   d2f = (f(:,2).*mob{3}(:,1) - f(:,1).*mob{3}(:,2) - ...
          2 .* df .* sum(mob{2}, 2)) ./ Lt;
end

%--------------------------------------------------------------------------

function varargout = mobilities(state, fluid)
   mu = fluid.properties(state);
   s  = fluid.saturation(state);
   [kr{1:nargout}] = fluid.relperm(s, state);

   %        \lambda_i in varargout{1}.
   % (d/ds) \lambda_i in varargout{2}.  Returned iff requested.
   %
   kr = cellfun(@(x)x(:,[1 end]), kr, 'UniformOutput', false);
   varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
                       'UniformOutput', false);
end
