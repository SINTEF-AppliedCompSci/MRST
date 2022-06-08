function fluid = makeVEFluid(Gt, rock, relperm_model, varargin)
%
% Construct a VE fluid with properties specific to a chosen model
%
% SYNOPSIS:
%   function fluid = makeVEFluid(Gt, rock, relperm_model, varargin)
%
% DESCRIPTION:
%
%
%
% PARAMETERS:
%   Gt            - Underlying top-surface grid with which the fluid object
%                   will be used.
%   rock          - Object holding the vertically-averaged rock properties of
%                   Gt (can be obtained from a normal rock structure using
%                   the 'averageRock' function in CO2lab).
%   relperm_model - Text string used to specify one of several possible
%                   models for computing upscaled permeabilities.  Options are:
%                   - 'simple' : sharp-interface model with linear relative
%                                permeabilities, no residual saturation, and
%                                vertically homogeneous rock
%                   - 'integrated' : sharp-interface model with linear
%                                    relative permeabilities.  Allows vertically
%                                    heterogeneous rock and impact of caprock
%                                    rugosity.
%                   - 'sharp interface' : sharp-interface model with linear relative
%                                         permeabilities and vertically
%                                         homogeneous rock.  Includes impact
%                                         of caprock rugosity.
%                   - 'linear cap.' : Linear capillary fringe model with
%                                     Brooks-Corey type relative
%                                     permeability and endpoint scaling.
%                   - 'S table' : capillary fringe model based on sampled
%                                 tables in the upscaled saturation parameter.
%                   - 'P-scaled table' : capillary fringe model based on sampled
%                                        tables in the upscaled capillary
%                                        pressure parameter.
%                   - 'P-K-scaled table' : capillary fringe model basd on
%                                          sampled tables in the upscaled
%                                          capillary pressure parameter, and
%                                          taking varations in permeability
%                                          into account through a Leverett
%                                          J-function relationship.
%                   A description of the different models can be found in the
%                   paper "Fully-Implicit Simulation of Vertical-Equilibrium
%                   Models with Hysteresis and Capillary Fringe" (Nilsen et
%                   al., Computational Geosciences 20, 2016).
%
%   varargin      - Optional arguments supplied as 'key'/value pairs
%                   ('pn'/nv).  These can be used to specify the dissolution
%                   model, subscale caprock rugosity and a range of other
%                   options.  See detailed documentation of available options
%                   in the function 'default_options()' below.
%
% OPTIONAL ARGUMENTS:
% A non-exhaustive overview of key optional arguments (refer to the internal
% function `default_options` to see the full range of options)
%
% Optional arguments related to the type of compressibility and viscosity
% model used
% 
%   fixedT      - If left empty, fluid properties will depend on pressure and
%                 temperature.  If set to a scalar, this will be considered
%                 the constant temperature of the simulated system, and fluid
%                 properties will only depend on pressure.
%  
%   co2_rho_ref - Reference density value for CO2 (used in black-oil formulation)
%   wat_rho_ref - Reference density value for water (used in black-oil formulation)
%   co2_rho_pvt - Compressibility model for CO2.  Possibilities are:
%                 * empty array ([]) - CO2 considered incompressible (uses
%                                      value for `co2_rho_ref`)
%                 * [cw, p_ref]      - constant compressibility.  `cw` is the
%                                      (scalar) compressibility, `p_ref` the
%                                      (scalar) reference pressure.
%                 * [pm, pM, tm, tM] - Interpolate from a sampled table that
%                                      covers the pressure interval [pm, pM] and
%                                      the temperature interval [tm, tM].
%                 The default option is number 3 on the above list.  A
%                 sampled table corresponding to the default values of
%                 `pm`/`pM` and `tm`/`tM` is provided with CO2lab.  Other
%                 tables can be generated with CoolProps.
%   wat_rho_pvt - Same as `co2_rho_pvt`, but for water/brine.
%   co2_mu_ref  - Reference viscosity for CO2
%   wat_mu_ref  - Reference viscosity for water
%   co2_mu_pvt  - Viscosity model for CO2.  Possibilities are:
%                 * empty array ([]) - CO2 viscosity considered constant
%                                      (uses value for `co2_mu_ref`)
%                 * [c, p_ref]       - Pressure-dependent viscosity with
%                                      constant coefficient (analog to
%                                      constant compressibility).  `c` is the
%                                      scalar coefficient, whereas `p_ref` is
%                                      the reference pressure.
%                 * [pm, pM, tm, tM] - Interpolate from a sampled table that
%                                      covers the pressure interval [pm, pM]
%                                      and the temperature interval [tm, tM].
%                 The default option is the first on the above list,
%                 i.e. constant viscostiy.
%   wat_mu_pvt  - Same as `co2_mu_pvt`, but for water/brine.
% 
% Optional arguments related to sampled property tables:
%   p_range / t_range - Each of these fields should be on the form [min, max], 
%                       descibing the pressure and temperature range for
%                       which the sampled property table should be generated.
%   pnum / tnum       - Number of (equidistant) samples for pressure and
%                       temperature when generating sampled property tables.
%
% Optional arguments related to residual saturation and dissolution 
% 
%  residual    - Two component vector, where first component represent
%                residual brine saturation and second component residual CO2
%                saturation.  Default is [0 0] (no residual saturation).
%  krmax       - Two component vector, representing fine-scale relative
%                permeabilities at end point saturation for brine (first
%                component) and CO2 (second component). Default will be set
%                to one minus the residual saturation of the opposite phase,
%                consistent with a fine-scale linear relative permeability
%                curve. This value is only relevant for the 'sharp
%                interface' relperm model (the other models are either using
%                linear relative permeabilities by default, or a capillary
%                fringe which requires the full CO2 relative permeability
%                curve to be provided through the optional 'kr3D' parameter).
%  dissolution - True or false, depending on whether or not to include CO2
%                dissolution into brine in the model.
%  dis_rate    - If dissolution is active, this option describes the
%                dissolution rate.  A zero value means "instantaneous"
%                dissolution, a positive value specifies a finite
%                rate. (default: 5e-11).
%  dis_max     - Maximum dissolution (default: 0.03)
%
% Optional arguments related to subscale rugosity of top surface
%   
%  surf_topo - What model to use for top surface rugosity.  Options are:
%              'smooth', 'sinus', 'inf_rough', and 'square'.  This option is
%              only used if the relperm model has been set to 'sharp
%              interface'.
%  top_trap  - Name of file containing the height of the top trap, either as
%              a single scalar value or as a value per cell.  Only used for
%              the relperm model 'sharp interface'.  Default is empty (no
%              subscale trapping).
%  surf_topo - Topography model used when computing the impact of caprock
%              rugosity.  Options are 'smooth', 'sinus', 'inf_rough', and
%              'square'.  Default is 'smooth'.
%     
% Optional arguments related to models for relative permeability and capillary
% pressure.  Relperm parameters are relevant for relperm-models 'S-table',
% 'P-scaled table', or 'P-K-scaled table'.
%
%  C               - scaling factor in Brooks-Corey type capillary pressure curve
%  alpha           - exponent used in Brooks-Corey type capillary pressure curve
%  beta            - exponent of Brooks-Corey type relperm curve
%  surface_tension - surface tension used in 'P-K-scaled table'
%  invPc3D         - inverse Pc function to use for computing capillary
%                    fringe.  If empt, a Brooks-Corey type curve will be
%                    constructed using 'C' and 'alpha' above.
%  kr3D            - CO2 relperm curve.  If empty, a Brooks-Corey type
%                    relperm curve with exponent 'beta' will be created.
%   
% Optional arguments related the rock matrix
%
%  pvMult_p_ref - Reference pressure for pore volume multiplier (default: 10 MPa)
%  pvMult_fac   - pore volume compressibility (default: 1e-5 / bar)
%  transMult    - modify transmissilbilties (such as due fo faults in the 3D grid)
%
% RETURNS:
%   fluid - struct containing the following functions (where X = 'W' [water]
%           and 'Gt' [gas])
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
%           * kr3D         - fine-scale relperm function
%           * invPc3D      - inverse fine-scale capillary pressure function
%
%           The following fields are optional, but may be returned by some
%           models:
%
%           * tranMultR(p) - mobility multiplier function
%           * transMult(p) - transmissibility multiplier function
%           * pvMult(p)    - pore volume multiplier function
%
% EXAMPLE:
%   The example script 'runStandardModel' (found under
%   'examples/publication_code/paper2') provides an example on how this
%   function is used.
%
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

   opt = merge_options(default_options(), varargin{:});
   fluid = []; % construct fluid from empty object

   %% Adding density and viscosity properties

   % Adding viscosity
   fluid = include_property(fluid, 'G', 'mu' , opt.co2_mu_ref,  opt.co2_mu_pvt ...
                            , opt.fixedT, opt.pnum, opt.tnum);
   fluid = include_property(fluid, 'W', 'mu' , opt.wat_mu_ref,  opt.wat_mu_pvt ...
                            , opt.fixedT, opt.pnum, opt.tnum);

   % Adding density
   fluid = include_property(fluid, 'G', 'rho', opt.co2_rho_ref, opt.co2_rho_pvt, ...
                            opt.fixedT, opt.pnum, opt.tnum);   
   fluid = include_property(fluid, 'W', 'rho', opt.wat_rho_ref, opt.wat_rho_pvt, ...
                            opt.fixedT, opt.pnum, opt.tnum);

   % Add density functions of the black-oil formulation type
   fluid = include_BO_form(fluid, 'G', opt.co2_rho_ref);
   fluid = include_BO_form(fluid, 'W', opt.wat_rho_ref);

   %% adding type-specific modifications
   fluid.surface_tension = opt.surface_tension;
   drho = fluid.rhoWS - fluid.rhoGS;
   C = opt.C * norm(gravity) * max(Gt.cells.H) * drho;
   switch relperm_model
     case 'simple' 
       fluid = setup_simple_fluid(fluid, Gt, opt.residual);
     case 'integrated' 
       fluid = setup_integrated_fluid(fluid, Gt, rock, opt.residual);
     case 'sharp interface'
       fluid = make_sharp_interface_fluid(fluid, Gt, opt.residual, opt.krmax, ...
                                          opt.top_trap, opt.surf_topo);
     case 'linear cap.'
       fluid = make_lin_cap_fluid(fluid, Gt, opt.residual);
     case 'S table'
       fluid = make_s_table_fluid(fluid, Gt, opt.residual, C, opt.alpha, ...
                                  opt.beta, opt.invPc3D, opt.kr3D);
     case 'P-scaled table'
       fluid = make_p_scaled_fluid(fluid, Gt, opt.residual, C, opt.alpha, ...
                                   opt.beta, opt.invPc3D, opt.kr3D);
     case 'P-K-scaled table'
       fluid = make_p_k_scaled_fluid(fluid, Gt, rock, opt.residual, C,  opt.alpha, ...
                                     opt.beta, opt.invPc3D, opt.kr3D);
     otherwise
       error([type, ': no such fluid case.']);
   end

   %% Adding dissolution-related modifications
   if opt.dissolution
      fluid.dis_rate = opt.dis_rate;
      fluid.dis_max  = opt.dis_max;
      fluid.rsSat    = @(pw, rs, flag, varargin) (pw*0+1) * fluid.dis_max;
   end

   %% Adding other modifications
   fluid.pvMultR = @(p) 1 + opt.pvMult_fac * (p - opt.pvMult_p_ref);

   % Transmissibility multipliers could be those supplied from the
   % topSurfaceGrid function to represent reduced transmissibilities in
   % partially overlapping faults.
   if ~isempty(opt.transMult)
      int_ix = ~any(Gt.faces.neighbors==0, 2);
      fluid.transMult = @(p) opt.transMult(int_ix);
   end
   
end

% ============================================================================

function opt = default_options()

   % Whether to include temperature as an argument in property functions
   opt.fixedT = []; % value of constant temperature field, or empty (if
                    % temperature should be an argument to the property
                    % functions.

   % Density of CO2 and brine
   opt.co2_rho_ref   =  760 * kilogram / meter^3; % Reference rho for CO2
   opt.wat_rho_ref   = 1100 * kilogram / meter^3; % Reference rho for brine

   p_range = [0.1, 400] * mega * Pascal; % CO2 default pressure range
   t_range = [  4, 250] + 274;           % CO2 default temperature range
   
   opt.pnum    = 800; % number of samples (if using EOS)
   opt.tnum    = 800; % number of samples (if using EOS)

   % The following options are used to specify whether CO2 and brine
   % densities should be considered constant, linear, or sampled from a
   % table.  If 'opt.co2_rho_pvt' is empty, CO2 will be considered as having
   % constant density (equal to the reference density provided in
   % 'opt.co2_rho_ref'). If it contains 2 values, [cw, p_ref], density will
   % be taken as a linear function of pressure, with cw being the
   % compressibility constant and 'p_ref' the pressure at which density
   % equals reference density.  If it contains four values, these will
   % represent the lower/upper bounds of pressure and temperature in the
   % sampled table to use.  The same explanation goes for water density,
   % as specified by 'opt.wat_rho_pvt'.
   opt.co2_rho_pvt = [p_range, t_range]; % empty, [cw, p_ref], or [pmin, pmax, tmin, tmax]
   opt.wat_rho_pvt = []; % empty, [cw, p_ref], or [pmin, pmax, tmin, tmax]

   % Viscosity of CO2 and brine

   % The explanation of the options 'co2_mu_pvt' and 'wat_mu_pvt' is
   % analoguous to that of 'co2_rho_pvt' and 'wat_rho_pvt' above.

   opt.co2_mu_ref = 6e-5 * Pascal * second;    % reference CO2 viscosity
   opt.wat_mu_ref = 8e-4 * Pascal * second;    % reference brine viscosity
   opt.co2_mu_pvt = []; % empty, [cw, p_ref], or [pmin, pmax, tmin, tmax]
   opt.wat_mu_pvt = []; % empty, [cw, p_ref], or [pmin, pmax, tmin, tmax]

   % Residual saturations [brine, co2]
   opt.residual = [0 0]; % default is no residual saturation for either phase
   opt.krmax = [-1 -1]; % endpoint relperms.  If negative (default), these
                        % will be set to one minus the residual saturation of
                        % the opposite phase.  This parameter is only
                        % relevant for the 'sharp interface' fluid model.

   % Dissolution of CO2 into brine
   opt.dissolution = false; % true or false
   opt.dis_rate    = 5e-11; % 0 means 'instantaneous'.  Otherwise, dissolution rate
   opt.dis_max     = 0.03;  % maximum dissolution

   % Caprock rugosity parameters (only used for the relperm model 'sharp
   % interface')
   opt.top_trap = [];
   opt.surf_topo = 'smooth'; % Choices are 'smooth', 'sinus', 'inf_rough',
                             % and 'square'.

   % Parameters used for relperms and capillary pressures based on table
   % lookup (i.e. 'S-table', 'P-scaled table' and 'P-K-scaled table)
   opt.C               = 0.4;   % scaling factor in Brooks-Corey type capillary pressure curve
   opt.alpha           = 0.5;   % Exponent used in Brooks-Corey type capillary pressure curve
   opt.beta            = 3;     % Exponent of Brooks-Corey type relperm curve
   opt.surface_tension = 30e-3; % Surface tension used in 'P-K-scaled table'
   opt.invPc3D         = [];    % Inverse Pc function to use for computing
                                % capillary fringe.  If empty, a Brooks-Corey
                                % type curve will be constructed using 'C', and
                                % 'alpha' above.
   opt.kr3D            = [];    % CO2 relperm curve.  If empty, Brooks-Corey
                                % type relperm curve with exponent 'beta'
                                % will be created.
   % Various parameters
   opt.pvMult_p_ref    = 100 * barsa;  % reference pressure for pore volume multiplier
   opt.pvMult_fac      = 1e-5 / barsa; % pore volume compressibility

   % Modify transmissibilities (such as due to faults in the 3D grid)
   opt.transMult = [];
   
end

% ----------------------------------------------------------------------------

function fluid = include_BO_form(fluid, shortname, ref_val)

   % Add reference value
   fluid.(['rho', shortname, 'S']) = ref_val;

   rhofun = fluid.(['rho', shortname]);

   if nargin(rhofun) > 1
      % density is considered function of P and T
      bfun = @(p, t, varargin) rhofun(p, t) ./ ref_val;
   else
      % density is considered a function of P only
      bfun = @(p, varargin) rhofun(p) ./ ref_val;
   end

   % Add fluid formation factor
   fluid.(['b', shortname]) = bfun;

end

% ----------------------------------------------------------------------------

function fluid = include_property(fluid, shortname, propname, prop_ref, prop_pvt, ...
                                  fixedT, pnum, tnum)
   if isempty(prop_pvt)
      % use constant property (based on reference property).  Whether it is a
      % function of P only, or of both P and T, is irrelevant here.
      fluid.([propname, shortname]) = as_function_of_p(prop_ref);

   elseif numel(prop_pvt) == 2
      % Linear compressibility
      cw    = prop_pvt(1);
      ref_p = prop_pvt(2);
      fluid.([propname, shortname]) = @(p, varargin) prop_ref * (1 + cw * (p - ref_p));

   else
      assert(isvector(prop_pvt) && numel(prop_pvt) == 4);

      fluid = addSampledFluidProperties(fluid, shortname, ...
                                        'pspan',  prop_pvt(1:2), ...
                                        'tspan',  prop_pvt(3:4), ...
                                        'pnum',   pnum, ...
                                        'tnum',   tnum, ...
                                        'props',  [strcmpi(propname, 'rho'), ...
                                                   strcmpi(propname, 'mu'),  ...
                                                   strcmpi(propname, 'h')], ...
                                        'fixedT', fixedT);
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

   fluid.krW = @(sw, varargin) sw;
   fluid.krG = @(sg, varargin) sg;
   fluid.kr3D = @(s, varargin) s;

end

% ----------------------------------------------------------------------------

function fluid = sharp_interface_cap_pressure(fluid, Gt)

   % When function is called, the following fields are needed
   require_fields(fluid, {'bW', 'bG', 'rhoWS', 'rhoGS'});

   % @@ Residual saturation not yet taken into account below!
   
   fluid.pcWG = @(sg, p, varargin)                       ...
                 norm(gravity) *                         ...
                 (fluid.rhoWS .* fluid.bW(p) -           ...
                  fluid.rhoGS .* fluid.bG(p)) .* (sg) .* ...
                  Gt.cells.H;

   fluid.invPc3D = @(p) 1 - (sign(p + eps) + 1) / 2;

end

% ----------------------------------------------------------------------------

function fluid = setup_simple_fluid(fluid, Gt, residual)

% Sharp interface; rock considered vertically uniform; no impact from caprock
% rugosity  

   fluid = linear_relperms(fluid);                    % 'krW'      , 'krG', 'kr3D'
   if residual(1) ~= 0 || residual(2) ~= 0
      warning('Simple model ignores nonzero residual saturation values.');
   end
   fluid.res_water = 0;
   fluid.res_gas = 0;
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

function fluid = make_sharp_interface_fluid(fluid, Gt, residual, krmax, top_trap, surf_topo)

% Sharp interface; rock considered vertically uniform; caprock rugosity
% influences relperm.  
   for i = 1:2
      if krmax(i) < 0
         krmax(i) = 1 - residual(3-i);
      end
   end
   fluid = addVERelperm(fluid       , Gt          , ...
                        'res_water' , residual(1) , ...
                        'res_gas'   , residual(2) , ...
                        'krw'       , krmax(1)    , ...
                        'krg'       , krmax(2)    , ...
                        'top_trap'  , top_trap    , ...
                        'surf_topo' , surf_topo);
end

% ----------------------------------------------------------------------------

function fluid = make_lin_cap_fluid(fluid, Gt, residual)

   % A model using a linear capillary fringe.  Rock considered vertically
   % uniform; no impact from caprock rugosity.

   % Local constants used:
   beta = 2;
   fac  = 0.2;
   g    = norm(gravity);
   drho = fluid.rhoWS - fluid.rhoGS;

   fluid = addVERelpermCapLinear(fluid, ...
                                 fac * g * max(Gt.cells.H) * drho , ...
                                 'res_water'   , residual(1)      , ...
                                 'res_gas'     , residual(2)      , ...
                                 'beta'        , beta             , ...
                                 'H'           , Gt.cells.H       , ...
                                 'kr_pressure' , true);
end

% ----------------------------------------------------------------------------

function fluid = make_s_table_fluid(fluid, Gt, residual, C, alpha, beta, ...
                                    invPc3D, kr3D)

   % Exact relationships

   % Local constants used
   drho    = fluid.rhoWS - fluid.rhoGS;
   samples = 1000;
   tabSw   = linspace(0, 1, 100)';
   tabW    = struct('S', 1 - tabSw, 'kr', tabSw, 'h', []);
   
   % Establish fine-scale functions
   if isempty(invPc3D)
      invPc3D = @(p) max((C ./ (p + C)).^(1 / alpha) , residual(1));
   end
   if isempty(kr3D)
      kr3D = @(s) s.^beta;
   end

   table_co2_1d = makeVEtables('invPc3D'    , invPc3D , ...
                               'is_kscaled' , false   , ...
                               'kr3D'       , kr3D    , ...
                               'drho'       , drho    , ...
                               'Gt'         , Gt      , ...
                               'samples'    , samples);

   fluid = addVERelperm1DTables(fluid , ...
        'res_water'   , residual(1)   , ...
        'res_gas'     , residual(2)   , ...
        'height'      , Gt.cells.H    , ...
        'table_co2'   , table_co2_1d  , ...
        'table_water' , tabW);
end

% ----------------------------------------------------------------------------

function fluid = make_p_scaled_fluid(fluid, Gt, residual, C, alpha, beta, ...
                                     invPc3D, kr3D)

   % Integral transformed from dz to dp

   % Local constants used
   drho    = fluid.rhoWS - fluid.rhoGS;
   samples = 2000;
   tabSw   = linspace(0, 1, 100)';
   tabW    = struct('S', 1 - tabSw, 'kr', tabSw, 'h', []);

   
   % Establish fine-scale functions
   if isempty(invPc3D)
      invPc3D = @(p) max((C ./ (p + C)).^(1 / alpha) , residual(1));
   end
   if isempty(kr3D)
      kr3D = @(s) s.^beta;
   end
   
   table_co2_1d = makeVEtables('invPc3D'    , invPc3D , ...
                               'is_kscaled' , false   , ...
                               'kr3D'       , kr3D    , ...
                               'drho'       , drho    , ...
                               'Gt'         , Gt      , ...
                               'samples'    , samples);

   fluid = addVERelperm1DTablesPressure(fluid , ...
       'res_water'   , residual(1)            , ...
       'res_gas'     , residual(2)            , ...
       'height'      , Gt.cells.H             , ...
       'table_co2'   , table_co2_1d           , ...
       'table_water' , tabW                   , ...
       'kr_pressure' , true);
 end

% ----------------------------------------------------------------------------

function fluid = make_p_k_scaled_fluid(fluid, Gt, rock, residual, C, alpha, beta, ...
                                       invPc3D, kr3D)

   % Integral transformed from dz to dp, and using Leverett's J-function

   % Local constants used
   kscale = sqrt(rock.poro ./ (rock.perm)) * fluid.surface_tension;
   drho    = fluid.rhoWS - fluid.rhoGS;
   samples = 2000;
   tabSw   = linspace(0, 1, 100)';
   tabW    = struct('S', 1 - tabSw, 'kr', tabSw, 'h', []);

   C = max(C ./ kscale);
   
   % Establish fine-scale functions
   if isempty(invPc3D)
      invPc3D = @(p) max((C ./ (p + C)).^(1 / alpha) , residual(1));
   end
   if isempty(kr3D)
      kr3D = @(s) s.^beta;
   end
   
   table_co2_1d = makeVEtables('invPc3D'    , invPc3D , ...
                               'is_kscaled' , true    , ...
                               'kr3D'       , kr3D    , ...
                               'drho'       , drho    , ...
                               'Gt'         , Gt      , ...
                               'samples'    , samples , ...
                               'kscale'     , kscale);

   fluid = addVERelperm1DTablesPressure(fluid , ...
       'res_water'   , residual(1)            , ...
       'res_gas'     , residual(2)            , ...
       'height'      , Gt.cells.H             , ...
       'table_co2'   , table_co2_1d           , ...
       'table_water' , tabW                   , ...
       'rock'        , rock                   , ...  % @@ full 3D rock?
       'kr_pressure' , true);

end

% ----------------------------------------------------------------------------
