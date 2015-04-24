function fluid = initSWOFFluid(varargin)
%Initialize incompressible two-phase fluid model (res. sat., an. rel-perm).
%
% SYNOPSIS:
%   fluid = initSWOFFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu     -- Phase viscosities in units of Pa*s.
%               - rho    -- Phase densities in units of kilogram/meter^3.
%               - table  -- deck.PROPS.SWOF (deck = readEclipseDeck(file)
%               - satnum -- deck.REGIONS.SATNUM
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


   opt = struct('table', [], 'mu', [], 'rho', [], 'satnum', -1);
   opt = merge_options(opt, varargin{:});

   opt.satnum = reshape(opt.satnum, [],1);
   %[krw, kro, Swco, krocw, Sorw, pcw] = swof(opt.table);
   [krw, kro, krocw, Swco, Swcr, Sowcr, Swmax, pcow] = swof(opt.table);
   prop = @(  varargin) properties(opt, varargin);


   if numel(krw) > 1 && numel(unique(opt.satnum))>1 %multiple regions
      assert(all(opt.satnum > 0))
      kr   = @(s,varargin) relperm_mult_regions(s, opt, krw, kro, Swco, ...
                                                krocw, Sowcr);
   else
      opt.satnum = max(opt.satnum(1), 1);
      kr  = @(s, varargin) relperm(s, opt, krw{opt.satnum}, kro{opt.satnum}, ...
                                  Swco(opt.satnum), krocw(opt.satnum), ...
                                  Sowcr(opt.satnum));
   end

   pc  = @(state) comp_capillary(state, pcow,  opt.satnum, Swco, Sowcr);

   fluid = struct('properties', prop             , ...
                  'saturation', @(x,varargin) x.s, ...
                  'relperm'   , kr, ...
                  'pc',@(x,varargin) pc(x));

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

function varargout = relperm(s, opt, krw, kro, Swco, krocw, Sowcr)
if(false)
   [s1, s2, den] = modified_saturations(s, [Swco, Sowcr]);

   [kw{1:nargout}] = krw(s1);
   [ko{1:nargout}] = kro(s2);
else
   [kw{1:nargout}] = krw(s(:,1));
   [ko{1:nargout}] = kro(1-s(:,1));

   den = 1;
end
   varargout{1}    = [ kw{1}, ko{1} ];

   if nargout > 1,
      null = zeros(size(kw{2}));
      varargout{2} = [kw{2}, null, null, ko{2}] ./ den;
   end
end

%--------------------------------------------------------------------------

function varargout = relperm_mult_regions(s, opt, krw, kro, Swco, krocw, Sorw)
if(size(s,1)~=numel(opt.satnum))
   % only allow use of non-matching satnum regions for estimating dt in
   % 'explicitTransport'
   st = dbstack;
   assert(any(cellfun(@(x) strcmp(x, 'estimate_dt'), {st.name})))
   %warning('saturation and satnum do not match we use region 1')
   opt.satnum=ones(size(s,1),1);
end
if(false)
   [s1, s2, den] = modified_saturations(s, [Swco(opt.satnum)'  ...
                                            Sorw(opt.satnum)']);

   [kw{1:nargout}] = evalMultipleRegions(krw, opt.satnum, s1);
   [ko{1:nargout}] = evalMultipleRegions(kro, opt.satnum, s2);
else
   [kw{1:nargout}] = evalMultipleRegions(krw, opt.satnum, s(:,1));
   [ko{1:nargout}] = evalMultipleRegions(kro, opt.satnum, 1-s(:,1));

   den = 1;
end
   varargout{1}    = [ kw{1}, ko{1} ];

   if nargout > 1,
      null = zeros(size(kw{2}));
      varargout{2} = [kw{2}, null, null, ko{2}] ./ den;
   end
end

function varargout = comp_capillary(state, table, satnum, Swco, Sorw)
if(false)
   s1 = modified_saturations(state.s, [Swco(satnum)' Sorw(satnum)']);
else
   s1 = state.s(:,1);
end

if numel(satnum)>1 % multiple regions
   [varargout{1:nargout}] = evalMultipleRegions(table, satnum, s1(:,1));
else
   [varargout{1:nargout}]= table{satnum}(s1); %[pc, dpc]
end
end


%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(s, sr)
   den = 1 - sum(sr, 2);
   s1  = (    s(:,1) - sr(:,1)) ./ den;  s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
   s2  = (1 - s(:,1) - sr(:,2)) ./ den;  s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
end
