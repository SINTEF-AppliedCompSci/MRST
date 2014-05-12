function fluid = initCoupledVEFluid(varargin)
% Initialize incompressible two-phase fluid model for vertical average
% calculation with both densities equal. This gives a simple realitic
% hysteres model with linear relperm functions
% 
%
% SYNOPSIS:
%   fluid = initCoupledVEFluid('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu  -- Phase viscosities in units of Pa*s.
%               - rho -- Phase densities in units of kilogram/meter^3.
%               - sr  -- Phase residual saturation
%               - height -- Height of all cells in the grid
%               - region3D -- numel(region) = g.cells.num, region3D = true for 3D cells.  
%               - n      -- Phase relative permeability exponents for 3D.
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
%   NB!  state has to have the fields s, s_max for this fluid to work 
%
% EXAMPLE:
%   fluid = initCoupledVEFluid('mu' , [   1,  10]*centi*poise     , ...
%                           'rho', [1014, 859]*kilogram/meter^3, ...
%                           'height'  , Gt.cells.H,...
%                           'sr', [0.2, 0.2], 'region3D', cells3D, 'n', [2 2]);
%
%   
%   s = linspace(0, 1, 1001).'; kr = fluid.relperm(state);
%   plot(s, kr), legend('kr_1', 'kr_2')
%
% SEE ALSO:
%   fluid_structure, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2012-10-04 14:50:34 +0200 (Thu, 04 Oct 2012) $
% $Revision: 10001 $

opt = struct('mu', [], 'rho', [],'sr', [], 'g', [],'region3D', [],'n', []);
opt = merge_options(opt, varargin{:});

%g_top = struct('H',opt.height);
g_top = opt.g;

n_mu = numel(opt.mu); n_rho = numel(opt.rho);n_sr=numel(opt.sr);
assert ((n_mu == 2) && (n_rho == 2) && (n_sr == 2));


prop = @(  varargin) properties(opt, varargin{:});
kr   = @(state) relperm(state,g_top,opt,varargin{:});
pc   = @(state) cap_press(state,g_top,opt, varargin{:});
sat2height =@(state)   saturation2Height(state,g_top,opt);
bndkr = @(s, smax, H) bndrelperm(s, smax, H, opt);
fluid = struct('properties', prop             , ...
               'saturation', @(x,varargin) x.s, ...%not used
               'relperm'   , @(s,state) kr(state),  ...
               'pc'        , @(state) pc(state),...
               'bndrelperm', @(s, smax, H) bndkr(s, smax, H), ...
               'sat2height', @(state)   sat2height(state));
end

%--------------------------------------------------------------------

function varargout = properties(opt, varargin)
   varargout{1}                 = opt.mu ;
   if nargout > 1, varargout{2} = opt.rho; end
   if nargout > 2, varargout{3} = opt.sr     ; end
end

%---------------------------------------------------------------------
function varargout = cap_press(state, g_top, opt, varargin)
   H = zeros(g_top.cells.num,1);
   H(~opt.region3D) = g_top.cells.H(g_top.cells.mapTopSurface(~opt.region3D));
   %(find(~opt.region3D)); %#ok
   % this trasformation has to be doen twice as long as
   % pc and relperm are in separate functions
   if nargout<2
      [h,h_max] = saturation2Height(state,g_top,opt);
   else
      [h,h_max,dh] = saturation2Height(state,g_top,opt);
   end
   
   ngrav=norm(gravity);
   %pc = zeros(numel(opt.region3D), 1);
   pc = ngrav*(bsxfun(@times,h,opt.rho(1))...
      +bsxfun(@times,H-h,opt.rho(2)));
   % pc(~opt.region3D) = pc_tmp(find(~opt.region3D)); %#ok
   pc(opt.region3D) = 0;
   varargout{1} = pc;
   if nargout==2
      % dpc = zeros(numel(opt.region3D), 1);
      dpc = ngrav*(opt.rho(1)-opt.rho(2))*dh;
      dpc(opt.region3D) = 0;
      varargout{2} = dpc;
   end
   if(nargout>2)
      error('3 output arguments not implemented');
   end
   
end
% NB: resdiual water saturation not properly implemented! kr(1-sw) = 1 now
function varargout = relperm(state, g_top, opt, varargin)
   %VE
   H = zeros(g_top.cells.num,1);
   H(~opt.region3D) = g_top.cells.H(g_top.cells.mapTopSurface(~opt.region3D));
   
   %(find(~opt.region3D)); %#ok
   if nargout<2
      [h,h_max] = saturation2Height(state,g_top,opt);
   else
      [h,h_max,dh] = saturation2Height(state,g_top,opt);
   end
   
   % kr = zeros(numel(opt.region3D), 2);
   %kr2D = bsxfun(@rdivide,[h, H-h],H);
   kr = bsxfun(@rdivide,[h, H-h],H);  %zeros(numel(opt.region3D), 2);
   %  kr(~opt.region3D, :) = kr2D(find(~opt.region3D),:); %#ok
   
   %3D
   s1 = state.s(opt.region3D,1); s2 = 1 - s1;
   kr(opt.region3D, :) = [s1 .^ opt.n(1), s2 .^ opt.n(2)];
   
   
   varargout{1} = kr;
   
   if nargout > 1,
      %dkr = zeros(numel(opt.region3D), 2);
      dkr = bsxfun(@rdivide,[dh,-dh],H);
      %dkr(~opt.region3D, :) = dkrtmp(find(~opt.region3D), :); %#ok
      dkr(opt.region3D, :) = [opt.n(1) .* s1 .^ (opt.n(1) - 1), ...
                             -opt.n(2) .* s2 .^ (opt.n(2) - 1)];
      varargout{2} = dkr;
   end
   if nargout > 2,
      varargout =zeros(numel(opt.region3D),2);
      %varargout =zeros(numel(H),2);
   end
   
end

function varargout = bndrelperm(s, s_max, H, opt, varargin)
   
   h = H.*(s-s_max.*(1-opt.sr(2)).*opt.sr(1))./(1-opt.sr(1)-opt.sr(2));
   h_max = s_max*H(1-opt.sr(2));
   
   kr = [h./H*(1-opt.sr(2)), (h_max-h)./H*(1-opt.sr(1)) + (H-h)*1./H];
   
   varargout{1} = kr;
   
   if nargout > 1,
      %does not account for residual saturation in CO2:
      dh=bsxfun(@rdivide,H,1-opt.sr(2)-opt.sr(1));
      %{
      if(h<h_max)
         dh=bsxfun(@rdivide,H,1-opt.sr(2)-opt.sr(1));
      else
         % this is not completely correct, it depends on sign ds
         % could be fixed by using sign of residual
         dh=bsxfun(@rdivide,H,1-opt.sr(2));
      end
      %}
      dkr = bsxfun(@rdivide,[dh,-dh], H);
      
      varargout{2} = dkr;
   end
   if nargout > 2,
      varargout =zeros(numel(opt.region),2);
      %varargout =zeros(numel(g_top.cells.H),2);
   end
end


function varargout = saturation2Height(state,g_top,opt)
   H = zeros(g_top.cells.num,1);
   H(~opt.region3D) = g_top.cells.H(g_top.cells.mapTopSurface(~opt.region3D));
   s_max = state.extSat(:,2);
   s     = state.s;
   assert(numel(s)==numel(H));
   h_max = bsxfun(@rdivide, s_max.*H, (1-opt.sr(2)) );
   h     = s.*H - bsxfun(@times,h_max,opt.sr(1));
   h     = bsxfun(@rdivide, h, (1-opt.sr(2)-opt.sr(1)) );
   %h = h(find(~opt.region3D)); %#ok
   %if nargout<3
   varargout{1}=h;
   varargout{2}=h_max;
   %end
   if nargout==3
      % if(h<h_max)
      %dh = zeros(
      dh=bsxfun(@rdivide,H,1-opt.sr(2)-opt.sr(1));
      % else
      % this is not completely correct, it depends on sign ds
      % could be fixed by using sign of residual
      % dh2(h>h_max)=bsxfun(@rdivide,H,1-opt.sr(2));
      
      %end
      %   dh = dh(find(~opt.region3D)); %#ok
      varargout{3}=dh;
   end
   if nargout>3
      error('wrong number of output')
   end
   % s*H=h*(1-sr(2))+(h_max -h)*sr(1)
   % s_max*H = h_max*(1-sr(2))
end
