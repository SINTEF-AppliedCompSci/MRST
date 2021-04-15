function fluid = initSimpleFluid(varargin)
%Initialize incompressible two-phase fluid model (analytic rel-perm).
%
% SYNOPSIS:
%   fluid = initSimpleFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters can be defined
%             with one value for each of the two fluid phases:
%               - mu  -- Phase viscosities in units of Pa*s. Defaults to 1
%                        1*centi*poise for each phase.
%               - rho -- Phase densities in units of kilogram/meter^3.
%                        Defaults to 1.
%               - n   -- Phase relative permeability exponents. 
%                        Defaults to 1.
%
% RETURNS:
%   fluid - Incompressible fluid data structure that is also compatible
%           with AD-based solvers.
%
% EXAMPLE:
%   fluid = initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
%                           'rho', [1014, 859]*kilogram/meter^3, ...
%                           'n'  , [   2,   2]);
%
%   s = linspace(0, 1, 1001).'; krW = fluid.krW(s); krO = fluid.krO(1 - s);
%   plot(s, [krW, krO]), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   `initSingleFluid`, `initSimpleFluidJfunc`, `incompTPFA`.

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

   opt = struct('mu' , [1, 1]*centi*poise, ... 
                'rho', [1, 1]*kilogram/meter^3, ...
                'n'  , [1, 1]);
   opt = merge_options(opt, varargin{:});

   assert(numel(opt.mu)  == 2, 'Viscosity (mu) must have exactly two entries.');
   assert(numel(opt.rho) == 2, 'Density (rho) must have exactly two entries.');
   assert(numel(opt.n)   == 2, 'Corey exponents (n) must have exactly two entries.');

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(s,varargin) relperm(s, opt, varargin{:});
   
   [bW , bO ] = recip_fvf();
   [krW, krO] = relperm_functions(opt);
   [muW, muO] = viscosity_functions(opt);

   fluid = struct('rhoWS', opt.rho(1), 'rhoOS', opt.rho(2), ...
                  'bW'   , bW        , 'bO'   , bO        , ...
                  'krW'  , krW       , 'krO'  , krO       , ...
                  'muW'  , muW       , 'muO'  , muO,        ...
                  'properties',         prop,               ...
                  'saturation',         @(x, varargin) x.s, ...
                  'relperm',            kr);
end

%--------------------------------------------------------------------------

function [bW, bO] = recip_fvf()
   [bW, bO] = deal(@(p) 1 + 0*p);  % Incompressible
end

%--------------------------------------------------------------------------

function [krW, krO] = relperm_functions(opt)
   krW = @(s) s .^ opt.n(1);
   krO = @(s) s .^ opt.n(2);
end

%--------------------------------------------------------------------------

function [muW, muO] = viscosity_functions(opt)
   muW = @(p) opt.mu(1) + 0*p;
   muO = @(p) opt.mu(2) + 0*p;
end

%--------------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = []     ; end
end

%--------------------------------------------------------------------------

function varargout = relperm(s, opt, varargin)
   s1 = s(:,1); s2 = 1 - s1;
   varargout{1} = [s1 .^ opt.n(1), s2 .^ opt.n(2)];

   if nargout > 1
      null = zeros([numel(s1), 1]);
      varargout{2} = [opt.n(1) .* s1 .^ (opt.n(1) - 1), ...
                      null, null                      , ...
                      opt.n(2) .* s2 .^ (opt.n(2) - 1)];
   end

   if nargout > 2
      a = opt.n .* (opt.n - 1);
      varargout{3} = [a(1) .* s1 .^ (opt.n(1) - 2), ...
                      a(2) .* s2 .^ (opt.n(2) - 2)];
   end
end
