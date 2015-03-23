function fluid = makeVEFluid(Gt, rock, relperm_model, varargin)
%
% Construct a VE fluid with properties specific to a chosen model
% 
% SYNOPSIS:
%   function fluid = makeADIFluid(type, opt, varargin)
%
% DESCRIPTION:
%
% @@ _complete_me_
%
% PARAMETERS:
%   Gt            - @@ 
%   relperm_model         - @@
%   varargin     - @@
%
% RETURNS:
%   fluid - struct containing the following functions (where X = 'W' [water]
%           and 'G' [gas])  
%           * rhoXS        - density of X at reference level (e.g. surface)
%           * bX(p), BX(p) - formation volume factors and their inverses
%           * muX(p)       - viscosity functions 
%           * krX(s)       - rel.perm for X
%           * rsSat(p)     - pressure-dependent max saturation value for
%                            dissolved gas 
%           * pcWG(sG, p)  - capillary pressure function
%           * dis_max      - maximum saturation value for dissolved gas
%           * dis_rate     - rate of dissolution of gas into water
%           * res_gas      - residual gas saturation
%           * res_water    - residual oil saturation
%           * kr3D         - @@
%           * invPc3D      - @@
%
%           The following fields are optional, but may be returned by some
%           models 
%           * tranMultR(p) - mobility multiplier function
%           * transMult(p) - transmissibility multiplier function
%           * pvMult(p)    - pore volume multiplier function
% EXAMPLE:
%
% SEE ALSO:
%
   
   opt = merge_options(default_options(), varargin{:});
   fluid = []; % construct fluid from empty object
   
   %% Adding density and viscosity properties
   
   % Adding viscosity
   fluid = include_property(fluid, 'G', 'mu' , opt.co2_mu_ref,  opt.co2_mu_pvt , opt.fixedT);
   fluid = include_property(fluid, 'W', 'mu' , opt.wat_mu_ref,  opt.wat_mu_pvt , opt.fixedT);

   % Adding density
   fluid = include_property(fluid, 'G', 'rho', opt.co2_rho_ref, opt.co2_rho_pvt, opt.fixedT);
   fluid = include_property(fluid, 'W', 'rho', opt.wat_rho_ref, opt.wat_rho_pvt, opt.fixedT);

   % Add density functions of the black-oil formulation type
   fluid = include_BO_form(fluid, 'G', opt.co2_rho_ref);
   fluid = include_BO_form(fluid, 'W', opt.wat_rho_ref);
   
   %% adding type-specific modifications
   switch relperm_model
     case 'simple' %@@ tested anywhere?
       fluid = setup_simple_fluid(fluid, Gt, opt.residual);
     case 'integrated' %@@ tested anywhere? 
       fluid = setup_integrated_fluid(fluid, Gt, rock, opt.residual);
     case 'sharp interface'
       fluid = make_sharp_interface_fluid(fluid, Gt, opt.residual);
     case 'linear cap.'
       fluid = make_lin_cap_fluid(fluid, Gt, opt.residual);     
     case 'S table'
       fluid = make_s_table_fluid(fluid, Gt, opt);
     case 'P-scaled table'
       fluid = make_p_scaled_fluid(fluid, Gt, opt);
     case 'P-K-scaled table'
       fluid = make_p_k_scaled_fluid(fluid, Gt, opt);
     otherwise
       error([type, ': no such fluid case.']);
   end
   
   %% Adding dissolution-related modifications
   if opt.dissolution
      f.dis_rate = opt.dis_rate;
      f.dis_max  = opt.dis_max;
      f.rsSat    = @(pw, rs, flag, varargin) (pw*0+1) * f.dis_max;
   end
            
   %% Adding other modifications
   f.pvMultR = @(p) 1 + opt.pvMult_fac * (p - opt.pvMult_p_ref);
   f.surface_tension = opt.surface_tension;
   
end

% ============================================================================

function opt = default_options()
   % Wheter to include temperature as an argument in property functions
   opt.fixedT = []; % value of constant temperature, or empty (if
                    % temperature should be an argument to the property
                    % functions. 

   % Density of CO2 and brine
   opt.co2_rho_ref   =  760 * kilogram / meter^3; % Reference rho for CO2
   opt.wat_rho_ref   = 1100 * kilogram / meter^3; % Reference rho for brine
   opt.co2_rho_pvt = []; % empty, function of (P) or of (P, T), or [pmin, pmax, tmin, tmax] 
   opt.wat_rho_pvt = []; % fct of (P, T); or name of sampled pv-table (if empty: constant rho)

   % Viscosity of CO2 and brine
   opt.co2_mu_ref = 6e-5 * Pascal * second;    % reference CO2 viscosity
   opt.wat_mu_ref = 8e-4 * Pascal * second;    % reference brine viscosity
   opt.co2_mu_pvt = []; % fct of (P, T); or name of sampled pv-table (if empty: constant mu)
   opt.wat_mu_pvt = []; % fct of (P, T); or name of sampled pv-table (if empty: constant mu)

   % Residual saturations [brine, co2]
   opt.residual = [0 0]; % default is no residual saturation for either phase
   
   % Dissolution of CO2 into brine
   opt.dissolution = false; % true or false
   opt.dis_rate    = 5e-11; % 0 means 'instantaneous'.  Otherwise, dissolution rate 
   opt.dis_max     = 0.03;  % maximum dissolution

   % Caprock rugosity parameters (only used for the relperm model 'sharp
   % interface')
   opt.top_trap = [];
   opt.surf_topo = 'smooth'; % Choices are 'smooth', 'sinus', 'inf_rough',
                             % and 'square'.
                                
   % Various parameters
   opt.pvMult_p_ref    = 100 * barsa;  % reference pressure for pore volume multiplier
   opt.pvMult_fac      = 1e-5 / barsa; % pore volume compressibility
   opt.surface_tension = 30e-3;        % Surface tension - used in some models

   
end

% ----------------------------------------------------------------------------

function fluid = include_BO_form(fluid, shortname, ref_val)
   
   % Add reference value
   fluid.(['rho', shortname, 'S']) = ref_val; 
   
   rhofun = fluid.(['rho', shortname]);
   
   if nargin(rhofun) > 1
      % density is considered function of P and T
      bfun = @(P, T) rhofun(P, T) ./ ref_val;
   else
      % density is considered a function of P only
      bfun = @(P) rhofun(P) ./ ref_val;
   end
   
   % Add fluid formation factor
   fluid.(['b', shortname]) = bfun;
      
end

% ----------------------------------------------------------------------------

function fluid = include_property(fluid, shortname, propname, prop_ref, prop_pvt, fixedT)

   if isempty(prop_pvt)
      % use constant property (based on reference property).  Whether it is a
      % function of P only, or of both P and T, is irrelevant here.
      fluid.([propname, shortname]) = as_function_of_p(prop_ref);
      
   elseif isa(prop_pvt, 'function_handle')
      
      % The exact function requested is already provided.  Just add it to
      % the fluid object.
      fluid.([propname, shortname]) = prop_pvt;
      
   else
      assert(isvector(prop_pvt) && numel(prop_pvt) == 4);

      fluid = addSampledFluidProperties(fluid, shortname, ...
                                        'pspan', prop_pvt(1:2), ...
                                        'tspan', prop_pvt(3:4), ...
                                        'props', [strcmpi(propname, 'rho'), ...
                                                  strcmpi(propname, 'mu'),  ...
                                                  strcmpi(propname, 'h')]);
   end
end

% ----------------------------------------------------------------------------

function fun = as_function_of_p(val)

   assert(isnumeric(val) && isscalar(val));
   fun = @(p, varargin) p * 0 + val;  
   
end

% ----------------------------------------------------------------------------

function require_fields(fluid, fields)

   % check that 'fluid' has all the required fields
   for f = fields
      assert(isfield(fluid, f{:}));
   end
end

% ----------------------------------------------------------------------------

function fluid = linear_relperms(fluid)

   fluid = setfield(fluid, 'krW' , @(sw, varargin) sg);
   fluid = setfield(fluid, 'krG' , @(sg, varargin) sg);
   % fluid = setfield(fluid, 'kr3D', @(s           )  s); @@ Used anywhere?
   
end

% ----------------------------------------------------------------------------

function fluid = residual_saturations(fluid, residual)

   fluid = setfield(fluid, 'res_water', residual(1));
   fluid = setfield(fluid, 'res_gas'  , residual(2));
   
end

% ----------------------------------------------------------------------------

function fluid = sharp_interface_cap_pressure(fluid, Gt)
   
   % When function is called, the following fields are needed
   require_fields(fluid, {'bW', 'bG', 'rhoWS', 'rhoGS'});
   
   fluid = setfield(fluid, 'pcWG', @(sg, p, varargin)                       ...
                                    norm(gravity) *                         ...
                                    (fluid.rhoWS .* fluid.bW(p) -           ...
                                     fluid.rhoGS .* fluid.bG(p)) .* (sg) .* ...
                                    Gt.cells.H); 
   
   fluid = setfield(fluid, 'invPc3D', @(p) 1 - (sign(p + eps) + 1) / 2);

end

% ----------------------------------------------------------------------------

function fluid = setup_simple_fluid(fluid, Gt, residual)
   
% Sharp interface; rock considered vertically uniform; no impact from caprock
% rugosity 
   
   fluid = linear_relperms(fluid);                    % 'krW'      , 'krG', 'kr3D'
   fluid = residual_saturations(fluid, residual);     % 'res_water', 'res_gas'
   fluid = sharp_interface_cap_pressure(fluid, Gt);    % 'pcWG'     , 'invPc3D'
   
end

% ----------------------------------------------------------------------------

function fluid = setup_integrated_fluid(fluid, Gt, rock, residual)
   
% Sharp interface; vertical variations in rock properties taken into account;
% caprock rugosity influences CO2 relperm
   
   fluid = addVERelpermIntegratedFluid(fluid         , ...
                                       'Gt'          , Gt          , ...
                                       'rock'        , rock        , ...
                                       'kr_pressure' , true        , ...
                                       'res_water'   , residual(1) , ...
                                       'res_gas'     , residual(2));
end

% ----------------------------------------------------------------------------

function fluid = make_sharp_interface_fluid(fluid, Gt, residual)

% Sharp interface; rock considered vertically uniform; caprock rugosity
% influences relperm
   
   fluid = addVERelperm(fluid       , Gt           , ...
                        'res_water' , residual(1)  , ...
                        'res_gas'   , residual(2)  , ...
                        'top_trap'  , opt.top_trap , ...
                        'surf_topo' , opt.surf_topo);
end

% ----------------------------------------------------------------------------

function fluid = make_lin_cap_fluid(fluid, G, residual)
   
   % A model using a linear capillary fringe.  Rock considered vertically
   % uniform; no impact from caprock rugosity.
      
   % Local constants used:
   beta = 2;
   fac  = 0.2;
   g    = norm(gravity);
   drho = fluid.rhoWS - fluid.rhoGS;
   
   fluid = addVERelpermCapLinear(fluid, ...
                     'res_water'   , residual(1)                      , ...
                     'res_gas'     , residual(2)                      , ...
                     'beta'        , 2                                , ...
                     'cap_scale'   , fac * g * max(Gt.cells.H) * drho , ...
                     'H'           , Gt.cells.H                       , ...
                     'kr_pressure' , true);
end

% ----------------------------------------------------------------------------

function fluid = make_s_table_fluid(fluid, G, opt)
   
   % Exact relationships
      
end

% ----------------------------------------------------------------------------

function fluid = make_p_scaled_fluid(fluid, G, opt)
   
   % Integral transformed from dz to dp
   
end

% ----------------------------------------------------------------------------

function fluid = make_p_k_scaled_fluid(fluid, G, opt)
   
   % Integral transformed from dz to dp, and using Leverett's J-function
   
end

% ----------------------------------------------------------------------------

