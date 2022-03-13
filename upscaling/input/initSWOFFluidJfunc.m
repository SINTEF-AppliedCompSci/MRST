function fluid = initSWOFFluidJfunc(varargin)
%Initialize incompressible two-phase fluid model with J scaling of pc.
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
%               - surf_tens -- surface tension in unit N/m
%               - jfunc  -- Whether or not pc given with J-scaling.
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
% EXAMPLE:
%
% SEE ALSO:
%   `fluid_structure`, `initSimpleFluid`, `solveIncompFlow`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('table', [], 'mu', [], 'rho', [], 'satnum', -1, ...
                'jfunc', false, 'rock', [], 'surf_tens', []);
   opt = merge_options(opt, varargin{:});

   assert(~opt.jfunc || (~isempty(opt.rock) && ~isempty(opt.surf_tens)))

   opt.satnum = reshape(opt.satnum, [],1);

   [krw, kro, krocw, Swco, Swcr, Sowcr, Swmax, pcow, pc_inv] = swof(opt.table);
   fwinv = compFracFlowInv(krw, kro, opt.mu, Swco, Sowcr);

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

   for i = 1:numel(opt.table)

      smax(i) = max(opt.table{i}(:,1));
      smin(i) = min(opt.table{i}(:,1));
   end

   smax = reshape(smax(opt.satnum), [],1);
   smin = reshape(smin(opt.satnum), [], 1);

   % nb: temporary use rock.perm(:,1), c++ code uses trace
   if opt.jfunc
       disp('Using rock.perm(:,1) for j-scaling');
       if size(opt.rock.perm, 2) ~= 3
       j_coeff = opt.surf_tens*sqrt(opt.rock.poro./opt.rock.perm(:,1));
       else
       j_coeff = opt.surf_tens*sqrt(opt.rock.poro./(sum(opt.rock.perm,2)/3));
       end
   else
      j_coeff = 1;
   end

   pc    = @(state) comp_capillary_jfunc(state, pcow, Swco, Sowcr, j_coeff, opt);
   pcinv = @(dp) pcinv_jfunc(dp,pc_inv,j_coeff, opt);
   fw_inv = @(f_w) f_inv(f_w, fwinv, opt);

   fluid = struct('properties', prop             , ...
                  'saturation', @(x,varargin) x.s, ...
                  'relperm'   , kr, ...
                  'pc',@(x,varargin) pc(x), ...
                  'pcinv', pcinv,...
                  'invpc', pcinv, ... % For upscaling project while refactoring.
                  'f_inv', fw_inv, ...
                  's_max', smax, ...
                  's_min', smin);
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
   %assert(any(cellfun(@(x) strcmp(x, 'estimate_dt'), {st.name})))
   if ~(any(cellfun(@(x) strcmp(x, 'estimate_dt'), {st.name})))
   warning('saturation and satnum do not match we use region 1')
   end
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

%--------------------------------------------------------------------------

function varargout = comp_capillary_jfunc(state, table, Swco, Sorw, j_coeff, opt)
if(false)
   s1 = modified_saturations(state.s, [Swco(opt.satnum)' Sorw(opt.satnum)']);
else
   s1 = state.s(:,1);
end

if numel(opt.satnum)>1 % multiple regions
   [varargout{1:nargout}] = evalMultipleRegions(table, opt.satnum, s1(:,1));
else
   [varargout{1:nargout}]= table{opt.satnum}(s1); %[pc, dpc]
end
for i = 1:numel(varargout)
   varargout{i} = varargout{i}.*j_coeff;
end
end
%--------------------------------------------------------------------------
function varargout = pcinv_jfunc(dp, table, j_coeff, opt)
j_func = dp./j_coeff;

if numel(opt.satnum)>1 % multiple regions
   [varargout{1:nargout}] = evalMultipleRegions(table, opt.satnum, j_func);
else
   [varargout{1:nargout}]= table{opt.satnum}(j_func); %[pc, dpc]
end
end
%--------------------------------------------------------------------------
function varargout = f_inv(f, table, opt)
if numel(opt.satnum)>1 % multiple regions
   [varargout{1:nargout}] = evalMultipleRegions(table, opt.satnum, f);
else
   [varargout{1:nargout}]= table{opt.satnum}(f);
end
end
%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(s, sr)
   den = 1 - sum(sr, 2);
   s1  = (    s(:,1) - sr(:,1)) ./ den;  s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
   s2  = (1 - s(:,1) - sr(:,2)) ./ den;  s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
end
