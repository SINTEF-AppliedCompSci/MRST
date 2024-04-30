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
%               - sr  -- Phase residual saturation [res_gas, res_water]
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
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('mu', [], 'rho', [],'sr', [0, 0],'height',[],'kwm',[1,1]);
opt = merge_options(opt, varargin{:});
g_top = struct('H',opt.height);

n_mu = numel(opt.mu); n_rho = numel(opt.rho);n_sr=numel(opt.sr);
assert ((n_mu == 2) && (n_rho == 2) && (n_sr == 2));
assert(numel(opt.height)>1);

prop = @(  varargin) properties(opt, varargin{:});
kr   = @(state) relperm(state,g_top,opt,varargin{:});
pc   = @(state) cap_press(state,g_top,opt, varargin{:});
% pc_inv=@(pc) cap_press_inv(pc,g_top,opt,varargin{:});
sat2height =@(state)   upscaledSat2height_wrapper(state,g_top,opt);
fluid = struct('properties', prop             , ...
               'saturation', @(x,varargin) x.s, ...%not used
               'relperm'   , @(s,state) kr(state),  ...
               'pc'        , @(state) pc(state),...
               'res_gas'   , opt.sr(1), ...
               'res_water' , opt.sr(2), ...
               'sat2height', @(state)   sat2height(state)); %,...
               %'invpc'     , @(pc)    pc_inv(pc),...
               % 's_max',1,...
               % 's_min',0);
               
end

%--------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = opt.sr ; end
end
%---------------------------------------------------------------------
% function s = cap_press_inv(pc, g_top, opt, varargin)
%    % this trasformation has to be done twice as long as
%    % pc and relperm are separate functions
%    %if nargout<2
%    %  [h,h_max] = upscaledSat2height_wrapper(state,g_top,opt);
%    %else
%    %  [h,h_max,dh] = upscaledSat2height_wrapper(state,g_top,opt);
%    %end
   
%    ngrav=norm(gravity);
   
%    %varargout{1} = ngrav*(bsxfun(@times,h,opt.rho(1))...
%    %   +bsxfun(@times,g_top.H-h,opt.rho(2)));
%    h=bsxfun(@rdivide,pc-ngrav*bsxfun(@times,g_top.H,opt.rho(2)),...
%        ngrav*(opt.rho(1)-opt.rho(2)));
%    % assuming no capillary fringe   
%    s=h./g_top.H;
%    s=max(s,0);
%    s=min(s,1);
% end
%---------------------------------------------------------------------
function varargout = cap_press(state, g_top, opt, varargin)
   % this trasformation has to be done twice as long as
   % pc and relperm are separate functions
   if nargout<2
     [h,h_max] = upscaledSat2height_wrapper(state,g_top,opt);
   else
     [h,h_max,dh] = upscaledSat2height_wrapper(state,g_top,opt);
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
     [h,h_max] = upscaledSat2height_wrapper(state,g_top,opt);
   else
     [h,h_max,dh] = upscaledSat2height_wrapper(state,g_top,opt);
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

function varargout = upscaledSat2height_wrapper(state, g_top, opt)
    sg = state.s;
    sgmax = state.extSat(:,2);
    
    % upscaledSat2height expect a grid, but as long as we do not do vertical
    % integration of heterogeneous pore volumes, we only need the grid
    % heights.  We create a mock grid object that only contains the cell
    % heights.
    Gt.cells.H = g_top.H;
    rg = opt.sr(1);
    rw = opt.sr(2);
    [h , hmax, dh] = upscaledSat2height(sg, sgmax, Gt, 'resSat', [rw, rg]);

    varargout = {h, hmax, dh};
end


