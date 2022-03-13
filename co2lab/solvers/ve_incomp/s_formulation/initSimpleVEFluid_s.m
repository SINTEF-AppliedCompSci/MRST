function fluid = initSimpleVEFluid_s(varargin)
% Initialize incompressible two-phase fluid model for vertical average
% calculation with both densities equal. This gives a simple realistic
% hysteresis model with linear relperm functions
% 
%
% SYNOPSIS:
%   fluid = initSimpleVEFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu  -- Phase viscosities in units of Pa*s.
%               - rho -- Phase densities in units of kilogram/meter^3.
%               - sr  -- Phase residual saturation
%               - height -- Height of all cells in the grid
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
%   NB!  state has to have the fields s, extSat for this fluid to work 
%
% EXAMPLE:
%   fluid = initSimpleVEFluid('mu' , [   1,  10]*centi*poise     , ...
%                           'rho', [1014, 859]*kilogram/meter^3, ...
%                           'height'  , Gt.cells.H,...
%                           'sr', [0.2, 0.2]);
%
%   
%   s = linspace(0, 1, 1001).'; kr = fluid.relperm(state);
%   plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   `fluid_structure`, `solveIncompFlow`.

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

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

opt = struct('mu', [], 'rho', [],'sr', [],'height',[],'kwm',[1,1]);
opt = merge_options(opt, varargin{:});
g_top = struct('H',opt.height);

n_mu = numel(opt.mu); n_rho = numel(opt.rho);n_sr=numel(opt.sr);
assert ((n_mu == 2) && (n_rho == 2) && (n_sr == 2));
assert(numel(opt.height)>1);

prop = @(  varargin) properties(opt, varargin{:});
kr   = @(state) relperm(state,g_top,opt,varargin{:});
pc   = @(state) cap_press(state,g_top,opt, varargin{:});
pc_inv=@(pc) cap_press_inv(pc,g_top,opt,varargin{:});
sat2height =@(state)   saturation2Height(state,g_top,opt);
fluid = struct('properties', prop             , ...
               'saturation', @(x,varargin) x.s, ...%not used
               'relperm'   , @(s,state) kr(state),  ...
               'pc'        , @(state) pc(state),...
               'invpc'     , @(pc)    pc_inv(pc),...
               'sat2height', @(state)   sat2height(state),...
               's_max',1,...
               's_min',0);
               
end

%--------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = opt.sr ; end
end
%---------------------------------------------------------------------
function s = cap_press_inv(pc, g_top, opt, varargin)
   % this trasformation has to be done twice as long as
   % pc and relperm are separate functions
   %if nargout<2
   %  [h,h_max] = saturation2Height(state,g_top,opt);
   %else
   %  [h,h_max,dh] = saturation2Height(state,g_top,opt);
   %end
   
   ngrav=norm(gravity);
   
   %varargout{1} = ngrav*(bsxfun(@times,h,opt.rho(1))...
   %   +bsxfun(@times,g_top.H-h,opt.rho(2)));
   h=bsxfun(@rdivide,pc-ngrav*bsxfun(@times,g_top.H,opt.rho(2)),...
       ngrav*(opt.rho(1)-opt.rho(2)));
   % assuming no capillary fringe   
   s=h./g_top.H;
   s=max(s,0);
   s=min(s,1);
end
%---------------------------------------------------------------------
function varargout = cap_press(state, g_top, opt, varargin)
   % this trasformation has to be done twice as long as
   % pc and relperm are separate functions
   if nargout<2
     [h,h_max] = saturation2Height(state,g_top,opt);
   else
     [h,h_max,dh] = saturation2Height(state,g_top,opt);
   end
   
   ngrav=norm(gravity);
   
   %varargout{1} = ngrav*(bsxfun(@times,h,opt.rho(1))...
   %   +bsxfun(@times,g_top.H-h,opt.rho(2)));
    varargout{1} = ngrav*(bsxfun(@times,h,opt.rho(1))...
       +bsxfun(@times,-h,opt.rho(2)));
   if nargout==2
     varargout{2} = ngrav*(opt.rho(1)-opt.rho(2))*dh;
   end
   if(nargout>2)
     error('3 output arguments not implemented');
   end
   
end
%---------------------------------------------------------------------
function varargout = relperm(state, g_top, opt, varargin)
   if nargout<2
     [h,h_max] = saturation2Height(state,g_top,opt);
   else
     [h,h_max,dh] = saturation2Height(state,g_top,opt);
   end
   
   %varargout{1} = bsxfun(@rdivide,[h*opt.kwm(1),(h_max-h)*opt.kwm(2)+(g_top.H-h_max)],g_top.H);
   varargout{1} = bsxfun(@rdivide,[h*opt.kwm(1),opt.kwm(2)*(g_top.H-h)],g_top.H);
   if nargout > 1,
      %varargout{2} = bsxfun(@rdivide,[dh*opt.kwm(1),-dh.*opt.kwm(2).*(h==h_max)-dh.*(h<h_max)],g_top.H);            
      varargout{2} = bsxfun(@rdivide,[dh*opt.kwm(1),-dh.*opt.kwm(2)],g_top.H);            
   end
   if nargout > 2,
     varargout =zeros(numel(g_top.H),2);
   end
end
function varargout = saturation2Height(state,g_top,opt)
   s_max = state.extSat(:,2);
   s     = state.s;
   assert(size(s,1)==numel(g_top.H));
   h_max = bsxfun(@rdivide, s_max.*g_top.H, (1-opt.sr(2)) );
   h     = s.*g_top.H - bsxfun(@times,h_max,opt.sr(1));
   h     = bsxfun(@rdivide, h, (1-opt.sr(2)-opt.sr(1)) );   
   %if nargout<3
%    assert(all(h>=0 | abs(h) < sqrt(eps)))
%    assert(all(h_max>=0))
%    assert(all(h<=g_top.H))
%    assert(all(h_max<=g_top.H))
   varargout{1}=h;
   varargout{2}=h_max;
   %end
   if nargout==3
     if(h<h_max)
       dh=bsxfun(@rdivide,g_top.H,1-opt.sr(2)-opt.sr(1));
     else
       % this is not completely correct, it depends on sign ds
       % could be fixed by using sign of residual
       dh=bsxfun(@rdivide,g_top.H,1-opt.sr(2));
     end
     varargout{3}=dh;
   end
   if nargout>3
     error('wrong number of output')
   end
   % s*H=h*(1-sr(2))+(h_max -h)*sr(1)
   % s_max*H = h_max*(1-sr(2))
end
