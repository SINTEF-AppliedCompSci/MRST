function fluid = initVEVerticalIntegratedFluid_s(varargin)
% Initialize incompressible two-phase fluid model for vertical average
% calculation, with bouth densities equal it wil be a simple realitic
% hysteres model with linear relperm functions
% 
%
% SYNOPSIS:
%   fluid = initVEVerticalIntegratedFluid_s('pn1', pv1, ...)
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
%   NB   state has to have the fields s,s_max for this fluid to work 
%
% EXAMPLE:
%   fluid = initVEVerticalIntegratedFluid('mu' , [   1,  10]*centi*poise     , ...
%                           'rho', [1014, 859]*kilogram/meter^3, ...
%                           'G'  , [],...
%                           'sr', [0.2, 0.2],...
%                            rockORG',[]);
%
%   
%   s = linspace(0, 1, 1001).'; kr = fluid.relperm(state);
%   plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   fluid_structure, solveIncompFlow.

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

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

   opt = struct('mu', [], 'rho', [],'sr', [],...
                'G',[],'int_poro',false,...
                'rockORG',[]);
   
   
   opt = merge_options(opt, varargin{:});
   assert(~isempty(opt.G));
   assert(~isempty(opt.rockORG));
   g_top = opt.G;
   if isfield(opt.G.cells,'H')
      g_top = opt.G;
      opt.H=g_top.cells.H;
      kr_H = integrateVertically(opt.rockORG.perm(:,1), g_top.cells.H, g_top);
      opt.perm2D=kr_H./g_top.cells.H;
      %kr_H=kr_H;./opt.perm2D(:,1);
      opt.kr_H=kr_H;
      pv_3D(g_top.columns.cells)=opt.rockORG.poro(g_top.columns.cells)...
      .*rldecode(g_top.cells.volumes,diff(g_top.cells.columnPos));
      opt.volumes= integrateVertically(pv_3D', inf, g_top).* ...
      g_top.cells.volumes;
   else
     %error('not implemented');
     %only need the columen structure and the dz structure
     if(isfield(opt.G,'columns'))
        G = opt.G;
        assert(isfield(G.cells,'columnPos'))
        assert(isfield(G.columns,'dz'))
     else 
        error('wrong grid')
     end
     kr_H = integrateVertically(opt.permORG(:,1), inf, g_top);
     opt.perm2D=kr_H./g_top.H;
     opt.kr_H=kr_H;
     opt.volumes=g_top.cells.volumes;
     tf=find(G.cells.faces(:,2)==6);
     assert(numel(tf)==G.cells.num)
     area=zeros(G.cells.num,1);
     cellno=rldecode(1:G.cells.num,diff(G.cells.facePos));
     area(cellno(tf))=G.faces.area(G.cells.faces(tf,1));
     opt.H=g_top.cells.volumes./area;
   end
   
   n_mu = numel(opt.mu); n_rho = numel(opt.rho);n_sr=numel(opt.sr);
   assert ((n_mu == 2) && (n_rho == 2) && (n_sr == 2));
   %assert(numel(opt.height)>1);

   prop = @(  varargin) properties(opt, varargin{:});
   kr   = @(state) relperm(state,g_top,opt,varargin{:});
   pc   = @(state) cap_press(state,g_top,opt, varargin{:});
   sat2height =@(state)   saturation2Height(state,opt);
   
   fluid = struct('properties', prop             , ...
                  'saturation', @(x,varargin) x.s, ...%not used
                  'relperm'   , @(s,state) kr(state),  ...
                  'pc'        , @(state) pc(state),...
                  'sat2height', @(state)   sat2height(state));
end

%--------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = []     ; end
end

%---------------------------------------------------------------------
function varargout = cap_press(state, g_top, opt, varargin)                %#ok
   % this trasformation has to be doen twise as long as
   % pc aand relperm i separate functions
   if nargout<2
      [h,h_max] = saturation2Height(state,opt);                            %#ok
     %[h,h] = saturation2Height(state,opt);
   else
     [h,h_max,dh] = saturation2Height(state,opt);                          %#ok
   end
   assert(all(h>=0));
   ng=norm(gravity);
   varargout{1} = ng*(bsxfun(@times,h,opt.rho(1))...
      +bsxfun(@times,opt.H-h,opt.rho(2)));
   if nargout==2
     varargout{2} =  ng*(opt.rho(1)-opt.rho(2))*dh;
   end
   if(nargout>2)
     error('3 output arguments not implemented');
   end
   
end

function varargout = relperm(state, g_top, opt, varargin)
   if opt.int_poro
     error('Int_poro: not implemented jet!')
%      if nargout<2
%        [h,h_max] = saturation2HeightIntporo(state,opt);
%      else
%        [h,h_max,dh] = saturation2HeightIntPoro(state,opt);
%      end
   else
     if nargout<2
        [h, h_max] = saturation2Height(state,opt);                         %#ok
     else
        [h, h_max, dh] = saturation2Height(state,opt);                     %#ok
     end
   end
   assert(all(h>=0));
   kr  = integrateVertically(opt.rockORG.perm(:,1), h, g_top);
   assert(all(kr>=0));
   varargout{1} = bsxfun(@rdivide,[kr, opt.kr_H-kr],opt.kr_H);
   if nargout > 1,
      varargout{2} = bsxfun(@rdivide,[dh,-dh],opt.H);            
   end
   if nargout > 2,
     varargout =zeros(numel(opt.H),2);
   end
end

function varargout = saturation2Height(state,opt)
   % this transformation is based on the simple transormation
   % s*H=h*(1-sr(2))+(h_max -h)*sr(1)
   % s_max*H = h_max*(1-sr(2))
   
   s_max=state.extSat(:,2);
   s = state.s;
   assert(all(s_max>=s));
   assert(numel(s)==numel(opt.H));
   h_max=bsxfun(@rdivide,s_max.*opt.H,(1-opt.sr(2)));
   h=s.*opt.H-bsxfun(@times,h_max,opt.sr(1));
   h=bsxfun(@rdivide,h,(1-opt.sr(2)-opt.sr(1)));
%   assert(all(h./opt.H>-eee) & all(h./opt.H<1+eee))
   h=max(h,0);
   h=min(h,opt.H);
   varargout{1}=h;
   varargout{2}=h_max;
   if nargout==3
     if(h<h_max)
       dh=bsxfun(@rdivide,opt.H,1-opt.sr(2)-opt.sr(1));
     else
       % this is not completely correct, it depends on sign ds
       % could be fixed by using sign of residual
       dh=bsxfun(@rdivide,opt.H,1-opt.sr(2));
     end
     varargout{3}=dh;
   end
   if nargout>3
     error('wrong number of output')
   end
end
function varargout = saturation2HeightIntPoro(state,opt)                   %#ok
   % this transformation is based on the simple tranformation
   % s*H=V(h)*(1-sr(2))+V(h_max) -V(h))*sr(1)
   % s_max*V(H) = V(h_max)*(1-sr(2))
   s_max=state.s_max;
   s = state.s;
   assert(numel(s)==numel(opt.H));
   Vh_max=bsxfun(@rdivide,s_max.*opt.volumes,(1-opt.sr(2)));
   V_h=s.*opt.volumes-bsxfun(@times,Vh_max,opt.sr(1));
   V_h=bsxfun(@rdivide,V_h,(1-opt.sr(2)-opt.sr(1)));
   h = opt.Vinv(V_h);
   h_max= opt.Vinv(Vh_max);
   if nargout<3
     varargout{1}=h;
     varargout{2}=h_max;
   end
   if nargout==3
     dh_ds=opt.volumes./opt.dV_dh(h);
     if(h<h_max)
       dh=bsxfun(@rdivide,dh_ds,1-opt.sr(2)-opt.sr(1));
     else
       % this is not completely correct, it depends on sign ds
       % could be fixed by using sign of residual
       dh=bsxfun(@rdivide,dh_ds,1-opt.sr(2));
     end
     varargout{3}=dh;
   end
   if nargout>3
     error('wrong number of output')
   end
end
