function fluid = initSingleFluid(varargin)
%Initialize incompressible single phase fluid model.
%
% SYNOPSIS:
%   fluid = initSingleFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined:
%               - mu    -- Fluid viscosity in units of Pa*s.
%               - rho   -- Fluid density in units of kilogram/meter^3.
%               - pahse -- Phase name.  Default value 'W'.  Supported
%               values are 'W', 'O', and 'G'.
%
% RETURNS:
%   fluid - Incompressible fluid data structure that is also compatible
%           with AD-based solvers.  Supports a single phase, identified by
%           the 'phase' parameter, only.  Requests for properties of other
%           phases return empty arrays/objects.
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

   opt = struct('mu', [], 'rho', [], 'phase', 'W');
   opt = merge_options(opt, varargin{:});

   assert (all([numel(opt.mu), numel(opt.rho)] == 1), ...
           'Function ''%s'' is supported for single phase only', mfilename);

   phase = upper(opt.phase(1));
   allph = {'W', 'O', 'G'};
   phIdx = strcmp(phase, allph);

   if ~any(phIdx)
      error('PhaseName:Unsupported', ...
           ['Phase name ''%s'' is unsupported.  Phase name ', ...
            'must be one of ''W'', ''O'', or ''G''.'], opt.phase);
   end

   phases = [ {phase}, allph(~phIdx) ];

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(s,varargin) relperm(s, opt, varargin{:});

   % 'P' is identified phase, 'M' is "Missing".
   [bP , bM ] = recip_fvf();
   [krP, krM] = relperm_functions();
   [muP, muM] = viscosity_functions(opt);

   fields = [ strcat('rho', phases, 'S'), ...
              strcat('b'  , phases),      ...
              strcat('kr' , phases),      ...
              strcat('mu' , phases),      ...
              { 'properties', 'saturation', 'relperm' }];

   values = [ { opt.rho, [] , []  }, ...
              { bP     , bM , bM  }, ...
              { krP    , krM, krM }, ...
              { muP    , muM, muM }, ...
              { prop, (@(x, varargin) x.s), kr } ];

   fluid = cell2struct(values, fields, 2);
end

%--------------------------------------------------------------------------

function [bP, bM] = recip_fvf()
   bP = @(p) 1 + 0*p;  % Incompressible
   bM = @(p) [];       % Missing
end

%--------------------------------------------------------------------------

function [krP, krM] = relperm_functions()
   krP = @(s) 1 + 0*s;  % Single phase => no relative permeability effects.
   krM = @(s) [];       % Missing
end

%--------------------------------------------------------------------------

function [muP, muM] = viscosity_functions(opt)
   muP = @(p) opt.mu + 0*p;
   muM = @(p) []; % Missing
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
