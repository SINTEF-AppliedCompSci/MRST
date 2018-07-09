function fluid = initSingleFluid(varargin)
%Initialize incompressible single phase fluid model.
%
% SYNOPSIS:
%   fluid = initSingleFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined:
%               - mu  -- Fluid viscosity in units of Pa*s.
%               - rho -- Fluid density in units of kilogram/meter^3.
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluid within the
%           reservoir model.
%
% EXAMPLE:
%   fluid = initSingleFluid('mu' ,    1*centi*poise     , ...
%                           'rho', 1014*kilogram/meter^3);
%
%   s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
%   plot(s, kr), legend('kr')
%
% SEE ALSO:
%   `fluid_structure`, `initSimpleFluid`, `solveIncompFlow`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


   opt = struct('mu', [], 'rho', []);
   opt = merge_options(opt, varargin{:});

   assert ((numel(opt.mu) == 1) && (numel(opt.rho) == 1));

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(s,varargin) relperm(s, opt, varargin{:});

   fluid = struct('properties', prop             , ...
                  'saturation', @(x,varargin) x.s, ...
                  'relperm'   , kr);
end

%--------------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1} = opt.mu;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = []     ; end
end

%--------------------------------------------------------------------------

function varargout = relperm(s, varargin)
   assert(size(s,2)==1);
   varargout{1} = ones(size(s));
   if nargout > 1, varargout{2} = zeros(size(s)); end
   if nargout > 2, varargout{3} = zeros(size(s)); end
end
