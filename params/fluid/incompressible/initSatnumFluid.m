function fluid = initSatnumFluid(T, varargin)
% Initialize incompressible two-phase fluid model with satnum and relperm
% evalution function.
%
%
% SYNOPSIS:
%   fluid = initSWOFFluid('pn1', pv1, ...)
%
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu     -- Phase viscosities in units of Pa*s.
%               - rho    -- Phase densities in units of kilogram/meter^3.
%
%             either: deck = readEclipseDeck(file)
%               - swof  -- table (deck.PROPS.SWOF)
%             or both of
%               - sgfn -- deck.PROPS.SGFN
%               - swfn -- deck.PROPS.SWFN
%
%               - satnum -- deck.REGIONS.SATNUM
%               - imbnum -- deck.REGIONS.IMBNUM
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
% EXAMPLE:
%   deck  = readEclipseDeck('TEST.DATA');
%   fluid = initSWOFFluid('rho',deck.PROPS.DENSITY([2,1]), ...
%           'mu', [deck.PROPS.PVTW(1), deck.PROPS.PVCDO(1)]*centi*poise,...
%           'table', deck.PROPS.SWOF, 'satnum', deck.REGIONS.SATNUM );
%
% SEE ALSO:
%   fluid_structure, initSimpleFluid, solveIncompFlow.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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



   opt = struct('mu', [], 'rho', [], ...
                'satnum', [], 'imbnum', []);

   opt = merge_options(opt, varargin{:});

   assert(isfield(T, 'swof')|| all([isfield(T, 'sgfn'), isfield(T, 'swfn')]));
   assert(~isempty(opt.satnum) && ~isempty(opt.imbnum))

   n_mu = numel(opt.mu); n_rho = numel(opt.rho);
   assert ((n_mu == 2) && (n_rho == 2));

   prop = @(  varargin) properties(opt, varargin);

   [kr, pc] = initSatnumRelPerm(T, opt.satnum, opt.imbnum);

   fluid = struct('properties', prop             , ...
                  'saturation', @(x, varargin) x.s, ...
                  'relperm'   , @(s, x, varargin) kr(s, x), ...
                  'pc', @(x,varargin) evalMultipleRegions(pc, opt.satnum, x.s(:,1)));
end

%--------------------------------------------------------------------------
function varargout = properties(opt, varargin)
      varargout{1} = opt.mu;
      if nargout > 1, varargout{2} = opt.rho; end
      if nargout > 2, varargout{3} = []     ; end
end

