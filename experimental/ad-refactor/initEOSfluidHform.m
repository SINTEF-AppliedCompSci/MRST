function fluid = initEOSfluidHform(Gt, EOSGAS, EOSWATER, varargin)
% Construct a fluid object that can be used with FullCompressibleCO2BrineModel
%
% SYNOPSIS:
%   function fluid = initEOSfluidHform(Gt, EOSGAS, EOSWATER, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt       - The top surface grid used (necessary for computing linear
%              relperms; not used if relperms are externally provided)
%   EOSGAS   - CO2 equation of state.  Must contain the functions 'rho(p,t)'
%              and 'mu(p,t)', _or_ constants 'rho' and 'mu'.  For fully
%              compressible model, higher derivatives is also needed.
%   EOSWATER - Water equation of state.  Must contain the functions rho(p,t)
%              and mu(p,t), _or_ constants 'rho' and 'mu'.
%   varargin - can be used to manually specify the relperm functions used
%              (otherwise, linear relperm functions will be provided).
%
% RETURNS:
%   fluid - fluid object with the fields:
%           fluid.gas.rho(p,t)   [gas density function]
%           fluid.gas.mu (p,t)   [gas viscosity function]
%           fluid.wat.rho(p,t)   [water densty function]
%           fluid.wat.mu (p,t)   [water viscosity function]
%           fluid.gas.kr(h)      [gas relperm function]
%           fluid.wat.kr(h)      [water relperm function, as a function of h_wat]
%
%           If used with a fully compressible model, the following higher
%           derivatives are needed (and must be provided in EOSGAS and/or
%           EOSWATER):
%
%           fluid.[gas|wat].beta(p,t)
%           fluid.[gas|wat].gamma(p,t)
%           fluid.[gas|wat].beta2(p,t)
%           fluid.[gas|wat].gamma2(p,t)
%           fluid.[gas|wat].chi(p,t)
%
%           If these fields are available, the following functions will
%           return non-empty values for the following functions:
%           
%           fluid.gas.h_integrals(p, t)
%           fluid.wat.h_integrals(p, t)
%           
% EXAMPLE:
%
% SEE ALSO:
%
   
   opt.krW                       = []; % function of 'h'. Empty value yields linear fun.
   opt.krG                       = []; % function of 'h'. Empty value yields linear fun.
   opt.constant_vertical_density = false;
   opt = merge_options(opt, varargin{:});   

   % Set density and viscosity functions
   [fluid.gas.rho, fluid.gas.mu] = deal(as_function(EOSGAS.rho),   as_function(EOSGAS.mu));
   [fluid.wat.rho, fluid.wat.mu] = deal(as_function(EOSWATER.rho), as_function(EOSWATER.mu));

   % Set relperm functions
   [fluid.gas.kr, fluid.wat.kr] = deal(opt.krG, opt.krW);
   
   if ~isa(fluid.gas.kr, 'function_handle') fluid.gas.kr = @(h_gas) h_gas ./ Gt.cells.H; end;
   if ~isa(fluid.wat.kr, 'function_handle') fluid.wat.kr = @(h_wat) h_wat ./ Gt.cells.H; end;
   
   % Set depth integrals for fully compressible model
   fluid.gas.h_integrals = setup_h_integral_fun(EOSGAS,   opt.constant_vertical_density);
   fluid.wat.h_integrals = setup_h_integral_fun(EOSWATER, opt.constant_vertical_density);
end

% ----------------------------------------------------------------------------
function res = setup_h_integral_fun(EOS, const_vertical_density)
% ----------------------------------------------------------------------------
% returns a function of 'slope' and 'Tgrad' that produces a set of functions
   EOS = add_extras(EOS);
   if complete_eos(EOS) && ~const_vertical_density
      % variable vertical density
      res =@(slope, tgrad) ...
           @(p, t) IetaAndINupEtaAndEta(p, t, EOS, ...
                                        tgrad/1000 * cos(slope), ...
                                        norm(gravity) * cos(slope));
   else
      % constant vertical density requested
      res = @(slope, tgrad) @(p, t) deal([], [], []);
   end
end

% ----------------------------------------------------------------------------
function [Ieta, INupEta, Eta] = IetaAndINupEtaAndEta(p, t, EOS, Gct, gct)
% ----------------------------------------------------------------------------
    EOS.compressible = 'full'; % required by the etaIntegrals function
    [Ieta, INupEta, ~, ~, Eta] = etaIntegrals(EOS, p , t, Gct, gct); 
end

% ----------------------------------------------------------------------------
function res = as_function(x) % ensure 'res' is a function of p and t
% ----------------------------------------------------------------------------
      if isa(x, 'function_handle') res = x; else res = @(p,t) 0 * p + x;  end
end
   
% ----------------------------------------------------------------------------
function res = complete_eos(EOS)
% ----------------------------------------------------------------------------
% Check if EOS has all the required functions to be useable for
% approximating vertical density profiles
   contains = @(name) isfield(EOS, name) && isa(EOS.(name), 'function_handle');
    
   % Return true if EOS contains all of the following functions:
   res = all(cellfun(contains, {'rho', 'beta', 'gamma', 'chi', 'beta2', 'gamma2'}));
end

% ----------------------------------------------------------------------------
function EOS = add_extras(EOS)
% ----------------------------------------------------------------------------
% add beta2, gamma2 and chi functions if they are not already there, and if
% the functions to construct them are avaialble.
   if isfield(EOS, 'rhoDPP') && ~isfield(EOS, 'beta2')
      EOS.beta2  = @(p,t) EOS.rhoDPP(p, t) ./ EOS.rho(p,t);
   end
   if isfield(EOS, 'rhoDTT') && ~isfield(EOS, 'gamma2')
      EOS.gamma2 = @(p,t) EOS.rhoDTT(p, t) ./ EOS.rho(p,t);
   end
   if isfield(EOS, 'rhoDPT') && ~isfield(EOS, 'chi')
      EOS.chi = @(p,t) EOS.rhoDPT(p, t) ./ EOS.rho(p, t);
   end
end
