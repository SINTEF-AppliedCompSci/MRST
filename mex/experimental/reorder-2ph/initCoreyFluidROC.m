function fluid = initCoreyFluidROC(varargin)
%Initialize incompressible two-phase fluid model compatible with reordering
%
% SYNOPSIS:
%   fluid = initCoreyFluidROC('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu  -- Phase viscosities in units of Pa*s.
%               - rho -- Phase densities in units of kilogram/meter^3.
%               - n   -- Phase relative permeability exponents.
%               - sr  -- Residual phase saturation.
%               - kwm -- Phase relative permeability at 'sr'.
%
%             Optionally, one may specify the following parameter
%               - reg -- Saturation regions, one value per cell in the
%                        model. NB! This option will only work if the
%                        fluid model is evaluated for each cell in the
%                        model at once.
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
% EXAMPLE:
%   fluid = initCoreyFluid('mu' , [   1,  10]*centi*poise     , ...
%                          'rho', [1014, 859]*kilogram/meter^3, ...
%                          'n'  , [   2,   2]                 , ...
%                          'sr' , [ 0.2, 0.2]                 , ...
%                          'kwm', [   1,   1]);
%
%   s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
%   plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   fluid_structure, initSimpleFluid, solveIncompFlow.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   opt = struct('mu', [], 'rho', [], 'n', [], 'sr', [], 'kwm', [], 'reg', []);
   opt = merge_options(opt, varargin{:});

   if ~isempty(opt.reg)
      m = 2*max(opt.reg);
      assert((numel(opt.n) >= m) & (numel(opt.sr) >= m) ...
         & (numel(opt.kwm) >= m));
   end

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(s,varargin) relperm(s, opt, varargin{:});

   fluid = struct('properties', prop             , ...
                  'saturation', @(x,varargin) x.s, ...
                  'relperm'   , kr               , ...
                  'param'     , opt);
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = []     ; end
end

%--------------------------------------------------------------------------

function varargout = relperm(s, opt, varargin)
   [s1, s2, den] = modified_saturations(s, opt);

   if isempty(opt.reg)
      n   = opt.n;
      kwm = opt.kwm;
   else
      assert(size(s,1)==numel(opt.reg));
      I = [2*opt.reg-1 2*opt.reg];
      n = opt.n(I);
      kwm = opt.kwm(I);
   end
   varargout{1}    = [ kwm(:,1) .* s1 .^ n(:,1), kwm(:,2) .* s2 .^ n(:,2)];

   if nargout > 1,
      null = zeros(size(s1));
      varargout{2} = [ kwm(:,1) .* n(:,1) .* s1 .^ (n(:,1) - 1), ...
                       null, null                       , ...
                       kwm(:,2) .* n(:,2) .* s2 .^ (n(:,2) - 1)] ./ den;
   end

   if nargout > 2,
      a = n .* (n - 1);
      varargout{3} = [ kwm(:,1) .* a(:,1) .* s1 .^ (n(:,1) - 2), ...
                       kwm(:,2) .* a(:,2) .* s2 .^ (n(:,2) - 2)] ./ den;
   end
end

%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(s, opt)
   if isempty(opt.reg)
      sr  = opt.sr;
   else
      assert(size(s,1)==numel(opt.reg));
      sr = opt.sr([2*opt.reg-1 2*opt.reg]);
   end

   den = 1 - sum(sr,2);
   s1  = (    s(:,1) - sr(:,1)) ./ den;  s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
   s2  = (1 - s(:,1) - sr(:,2)) ./ den;  s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
end
