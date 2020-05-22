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
%   fluid - Incompressible fluid data structure that is also compatible
%           with AD-based solvers.  Supports a 'water' phase only.
%           Requests for 'oil' properties return empty arrays/objects.
%
% EXAMPLE:
%   fluid = initSingleFluid('mu' ,    1*centi*poise     , ...
%                           'rho', 1014*kilogram/meter^3);
%
%   s = linspace(0, 1, 1001).'; kr = fluid.krW(s);
%   plot(s, kr), legend('kr')
%
% SEE ALSO:
%   `initSimpleFluid`, `incompTPFA`.

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

   opt = struct('mu', [], 'rho', []);
   opt = merge_options(opt, varargin{:});

   assert (all([numel(opt.mu), numel(opt.rho)] == 1), ...
           'Function ''%s'' is supported for single phase only', mfilename);

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(s,varargin) relperm(s, opt, varargin{:});

   [bW , bO ] = recip_fvf();
   [krW, krO] = relperm_functions();
   [muW, muO] = viscosity_functions(opt);

   fluid = struct('rhoWS', opt.rho, 'rhoOS', [] ,       ...
                  'bW'   , bW     , 'bO'   , bO ,       ...
                  'krW'  , krW    , 'krO'  , krO,       ...
                  'muW'  , muW    , 'muO'  , muO,       ...
                  'properties',     prop,               ...
                  'saturation',     @(x, varargin) x.s, ...
                  'relperm',        kr);
end

%--------------------------------------------------------------------------

function [bW, bO] = recip_fvf()
   bW = @(p) 1 + 0*p;  % Incompressible
   bO = @(p) [];       % Missing
end

%--------------------------------------------------------------------------

function [krW, krO] = relperm_functions()
   krW = @(s) 1 + 0*s;  % Single phase => no relative permeability effects.
   krO = @(s) [];       % Missing
end

%--------------------------------------------------------------------------

function [muW, muO] = viscosity_functions(opt)
   muW = @(p) opt.mu + 0*p;
   muO = @(p) []; % Missing
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
