function fluid = initSimpleVEFluidSFormSara1D(varargin)
% Initialize incompressible two-phase fluid model for vertical average
% calculation with both densities equal. This gives a simple realitic
% hysteres model with linear relperm functions
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
%   fluid_structure, solveIncompFlow.

%{
#COPYRIGHT#
%}

% $Date: 2012-01-30 11:41:03 +0100 (Mon, 30 Jan 2012) $
% $Revision: 9020 $

opt = struct('mu', [], 'rho', [],'sr', [],'height',[],'table_co2',[],'table_water',[]);
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
   % this trasformation has to be doen twise as long as
   % pc aand relperm i separate functions
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
   
   if(false)
    h_s=h./g_top.H;   
    s = interpTable(opt.table_co2.h,opt.table_co2.S,h_s);
   else
    sH = interpTable(opt.table_co2.hH,opt.table_co2.SH,h);
    s  = sH./g_top.H;
   end
end
%---------------------------------------------------------------------
function varargout = cap_press(state, g_top, opt, varargin)
   % this trasformation has to be doen twise as long as
   % pc aand relperm i separate functions
   if(false)
    h_s=interpTable(opt.table_co2.S,opt.table_co2.h,state.s(:,1));
    dh_s=dinterpTable(opt.table_co2.S,opt.table_co2.h,state.s(:,1));
    h=h_s.*g_top.H;
    dh=dh_s.*g_top.H;
   else
     SH= state.s(:,1).*g_top.H; 
     h=interpTable(opt.table_co2.SH,opt.table_co2.hH,SH);
     dh=dinterpTable(opt.table_co2.SH,opt.table_co2.hH,SH).*g_top.H;  
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
function varargout = relperm(state,g_top, opt, varargin)
   kr=nan(size(state.s,1),2);
   dkr=zeros(size(state.s,1),4);
   kr(:,1)=interpTable(opt.table_co2.S,opt.table_co2.kr,state.s(:,1));
   dkr(:,1)=dinterpTable(opt.table_co2.S,opt.table_co2.kr,state.s(:,1));
   kr(:,2)=interpTable(opt.table_water.S,opt.table_water.kr,1-state.s(:,1));
   dkr(:,4)=dinterpTable(opt.table_water.S,opt.table_water.kr,1-state.s(:,1)); 
   varargout{1} = kr;
   varargout{2} = dkr;
end

