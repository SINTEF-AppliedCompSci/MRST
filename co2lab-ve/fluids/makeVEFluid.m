function fluid = makeVEFluid(Gt, rock, relperm_model, varargin)
%
% Construct a VE fluid with properties specific to a chosen model
%
% SYNOPSIS:
%   function fluid = makeVEFluid(Gt, rock, relperm_model, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt            - Underlying top-surface (vertically averaged, 2D) grid with 
%                   which the fluid object will be used.
%   rock          - Object holding the vertically-averaged (i.e. 2D) rock
%                   properties associated with Gt. Can be obtained from the
%                   corresponding 3D rock structure using the `averageRock` 
%                   function in co2lab-common.
%   relperm_model - Text string used to specify one of several possible
%                   models for computing upscaled permeabilities.  Options
%                   are:
%                   - 'sharp_interface_simple': sharp interface model with
%                                               vertically averaged
%                                               properties, yielding linear
%                                               upscaled relperm curves.
%                   - 'sharp_interface_integrated': sharp interface model
%                                                   with integration of vertical
%                                                   heterogeneities, yielding
%                                                   nonlinear relperm curves.
%                   - 'P-scaled table' : capillary fringe model based on sampled
%                                        tables in the upscaled capillary
%                                        pressure parameter.
%                   - 'P-K-scaled table' : capillary fringe model based on
%                                          sampled tables in the upscaled
%                                          capillary pressure parameter, and
%                                          taking varations in permeability
%                                          into account through a Leverett
%                                          J-function relationship.
%                   - 'S table' : capillary fringe model based on sampled
%                                 tables in the upscaled saturation
%                                 parameter.  Only valid for aquifers with
%                                 constant height, and is provided here
%                                 mainly for reference purposes.
%                   The three capillary fringe models come with the following 
%                   limitations:
%                   - vertical rock heterogeneity is assumed
%                   - they do not combine with upscaled caprock rugosity (see
%                   'rugosity' parameter further down).
%                   If you choose to use a capillary fringe model, you also
%                   need to specify the fine-scale permeability ('kr3D') and
%                   inverse capillary ('invPc3D') functions.  See further
%                   down for details on how to specify these parameters.
%                   
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
%   reservoirT  - Reservoir temperature.  Used in the computation of fluid 
%                 properties when sampled tables are used (cf. `co2_rho_pvt`,
%                 `co2_mu_pvt`, `wat_rho_pvt` and `wat_mu_pvt` below).  
%                 If sampled tables are not used for fluid properties, then 
%                 the value of `reservoirT` will have no effect.  
%                 This option can be set to a scalar value (one single
%                 temperature for the whole reservoir), or a vector of values
%                 (one value per cell).  The latter is useful when there is a
%                 range of temperatures in the reservoir, e.g. a large
%                 sloping reservoir in the presence of a geothermal gradient.
%                 Default value is 30 degrees C (303.15 degrees Kelvin).
%  
%   co2_rho_ref - Reference density value for CO2 (used in black-oil formulation)
%   wat_rho_ref - Reference density value for water (used in black-oil formulation)
%   co2_rho_pvt - Compressibility model for CO2.  Possibilities are:
%                 (1) empty array ([]) - CO2 considered incompressible (uses
%                                        value for `co2_rho_ref`)
%                 (2) [cw, p_ref]      - constant compressibility.  `cw` is the
%                                        (scalar) compressibility, `p_ref` the
%                                        (scalar) reference pressure.
%                 (3) [pm, pM, tm, tM] - Interpolate from a sampled table that
%                                        covers the pressure interval [pm, pM] and
%                                        the temperature interval [tm, tM].
%                 The default option is (3) above.  A sampled table
%                 corresponding to the default values of `pm`/`pM` and
%                 `tm`/`tM`  is provided with CO2lab.  Other tables can be
%                 generated on the fly (requires `CoolProp` installed).
%                 If you use a sampled table, make sure that it covers all
%                 the pressures that might be reached, and all temperatures
%                 encountered in `reservoirT`.
% 
%   wat_rho_pvt - Same as `co2_rho_pvt`, but for water/brine.
%   co2_mu_ref  - Reference viscosity for CO2
%   wat_mu_ref  - Reference viscosity for water
%   co2_mu_pvt  - Viscosity model for CO2.  Possibilities are:
%                 (1) empty array ([]) - CO2 viscosity considered constant
%                                        (uses value for `co2_mu_ref`)
%                 (2) [c, p_ref]       - Pressure-dependent viscosity with
%                                        constant coefficient (analog to
%                                        constant compressibility).  `c` is the
%                                        scalar coefficient, whereas `p_ref` is
%                                        the reference pressure.
%                 (3) [pm, pM, tm, tM] - Interpolate from a sampled table that
%                                        covers the pressure interval [pm, pM]
%                                        and the temperature interval [tm, tM].
%                 The default option is the first on the above list,
%                 i.e. constant viscosity. If you use a sampled table, make
%                 sure that it covers all the pressures that might be reached,
%                 and all temperatures encountered in `reservoirT`.
%
%   wat_mu_pvt  - Same as `co2_mu_pvt`, but for water/brine.
% 
% Optional arguments related to sampled property tables:
%   pnum / tnum       - Number of (equidistant) samples for pressure and
%                       temperature when generating sampled property tables.
%                       (Default value is 800 for both).
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
%                curve. This value is only relevant for the sharp interface
%                relperm models, i.e. those with no capillary fringe.
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
%  rugosity  - scalar or field with one value per cell in the topsurface
%              grid, representing the upscaled structural trapping potential
%              within the cell.  It is given as a length, which expresses
%              the thickness of the "accretion layer" representing the
%              effect of the subscale traps.  Rugosity can be used with the
%              sharp interface relperm models, but does not combine with the
%              capillary fringe models (in which case it is ignored).
%  
% If a capillary fringe model is chosen, you also need to specify the
% fine-scale relative permeabilities and capillary pressure.  This is done
% with the following arguments:
% kr3D - CO2 relperm curve.  Either specified as a function handle (function
%        of saturation), or as scalar exponent of a Brooks-Corey type relperm
%        curve, i.e. the 'beta' in a function like kr = max(s-src, 0).^beta,
%        where 'src' is the residual CO2 saturation.  Default is 3 (cubic relperm).
% invPc3D - Inverse capillary pressure function used for computing the
%           capillary fringe saturation profile.  Either specified as a
%           function handle (saturation as a function of capillary pressure),
%           or as the scalar pair [C, alpha], which respectively represent
%           the scaling factor and the exponent for an inverse Brooks-Corey 
%           capillary pressure curve, 
%           i.e. sw = max( (C* ./ (p + C*)).^(1 / alpha), srw), where 'srw' is
%           the residual water saturation, and C* equals C scaled by the
%           factor max(Gt.cells.H) * norm(gravity) * (wat_rho_ref - CO2_rho_ref).
% 
%           Note that if the 'P-K-Scaled Table' relperm model is chosen, then
%           the function specified by 'invPc3D' represents the inverse of the
%           Everett J-function, which will be properly scaled to produce the
%           inverse capillary pressure based on the rock properties of each
%           vertical grid column.  
%           The default value of this parameter is [0.4, 0.5].
% surface_tension - This parameter is used when the 'P-K-Scaled Table'
%                   relperm model is chosen.  It represents the fluid surface
%                   tension, multiplied by the cosine of the contact angle,
%                   and is used along with rock properties to compute the
%                   multiplicative factor that converts the Everett
%                   J-function into a (rock property dependent) capillary
%                   pressure.  If any other relperm model is chosen, this
%                   variable is ignored.  Default value is 0.03.
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
%   The example script 'exampleVE' provides an example on how this
%   function is used.
%
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

   opt = merge_options(default_options(), varargin{:});
   if isempty(opt.invPc3D)
       opt.invPc3D = [0.4, 0.5];    
   end

   fluid = []; % construct fluid from empty object

   %% Adding density and viscosity properties

   % Adding viscosity
   fluid = include_property(fluid, 'G', 'mu' , opt.co2_mu_ref,  opt.co2_mu_pvt ...
                            , opt.reservoirT, opt.pnum, opt.tnum);
   fluid = include_property(fluid, 'W', 'mu' , opt.wat_mu_ref,  opt.wat_mu_pvt ...
                            , opt.reservoirT, opt.pnum, opt.tnum);

   % Adding density
   fluid = include_property(fluid, 'G', 'rho', opt.co2_rho_ref, opt.co2_rho_pvt, ...
                            opt.reservoirT, opt.pnum, opt.tnum);   
   fluid = include_property(fluid, 'W', 'rho', opt.wat_rho_ref, opt.wat_rho_pvt, ...
                            opt.reservoirT, opt.pnum, opt.tnum);

   % Add density functions of the black-oil formulation type
   fluid = include_BO_form(fluid, 'G', opt.co2_rho_ref);
   fluid = include_BO_form(fluid, 'W', opt.wat_rho_ref);
   
   %% Adding residual saturations
   fluid.res_water = opt.residual(1);
   fluid.res_gas = opt.residual(2);

   %% adding type-specific modifications
   fluid.surface_tension = opt.surface_tension;

   krw = ifelse(opt.krmax(1) < 0, 1 - opt.residual(2), opt.krmax(1));
   krg = ifelse(opt.krmax(2) < 0, 1 - opt.residual(1), opt.krmax(2));

   fluid.relperm_model = relperm_model;
   switch relperm_model
     case 'sharp_interface_simple'
       fluid = addVERelpermSharpInterface(fluid, Gt, rock, ...
                                          'dh', opt.rugosity, ...
                                          'krw', krw, ...
                                          'krg', krg, ...
                                          'type', 'simple');
     case 'sharp_interface_integrated'
       fluid = addVERelpermSharpInterface(fluid, Gt, rock, ...
                                          'dh', opt.rugosity, ...
                                          'krw', krw, ...
                                          'krg', krg, ...
                                          'type', 'integrated');
     case {'S table', 'P-scaled table', 'P-K-scaled table'}
       fun3D = setup_fine_scale_functions(opt.invPc3D, opt.kr3D, fluid, Gt);
       fluid = addVERelpermCapillaryFringe(fluid, Gt, rock, fun3D.invPc3D, ...
                                           fun3D.kr3D, 'type', relperm_model);
       if isfield(fun3D, 'pc3D')
           % the fine-scale capillary function (of sw) is given.  Store it
           % for reference.
           fluid.pc3D = fun3D.pc3D;
       end
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

% ----------------------------------------------------------------------------
function fun3D = setup_fine_scale_functions(invPc3D, kr3D, fluid, Gt)

    
    % setup fine-scale CO2 relperm function
    if isempty(kr3D)
        kr3D = 3; % default value (but cannot be set directly in default_options)
    end
    if isnumeric(kr3D)
        assert(isscalar(kr3D)); % should be a single number
                                % define a simple Corey
        beta = kr3D; % a scalar exponent was provided
        kr3D = @(s) max(s - fluid.res_gas, 0) .^ beta;
    end
    assert(isa(kr3D, 'function_handle'));
    if kr3D(fluid.res_gas/2) > 0
        warning(['Provided fine-scale relperm function for CO2 inconsistent ' ...
                 'with value provided for residual saturation.']);
    end
    
    fun3D.kr3D = kr3D;
    % setup fine-scale inverse capillary pressure function
    if isnumeric(invPc3D)
        assert(numel(invPc3D) == 2); % [C, alpha]
        
        drho = fluid.rhoWS - fluid.rhoGS;
        Hmax = max(Gt.cells.H);
        [C, a] = deal(invPc3D(1), invPc3D(2));
        C = C * norm(gravity) * Hmax * drho;
        
        invPc3D = @(p) max( (C ./ (p + C)).^(1 / a), fluid.res_water);
        
        % we store the corresponding fine-scale capillary function too, for
        % reference (though it is not used in the VE simulation)
        fun3D.pc3D = @(sw) C .* (1./max(sw, fluid.res_water).^a - 1);
        
    end
    assert(isa(invPc3D, 'function_handle'));
    fun3D.invPc3D = invPc3D;
end

% ============================================================================

function opt = default_options()

   % Whether to include temperature as an argument in property functions
   opt.reservoirT = 273.15 + 30; % (constant) reservoir temperature.  Can be specified
                             % as a scalar giving the same temperature for
                             % the whole reservoir, or as a vector with one
                             % temperature per gridcell (useful to model
                             % varying temperature in a sloping aquifer due
                             % to the thermal gradient.

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
   opt.rugosity = 0;
   
   % Parameters used for relperms and capillary pressures based on table
   opt.surface_tension = 30e-3; % Surface tension used in 'P-K-scaled table'
   opt.kr3D = [];                % fine-scale relative permeability (scalar
                                 % or function handle).  Since types of the
                                 % two options differ, it cannot be set as a
                                 % default here.
   opt.invPc3D = []; % Inverse fine-scale capillary pressure function.  This
                     % is specified either directly as a function handle, or
                     % as the scalar pair [C, alpha], which represent
                     % the scaling factor and the exponent for an inverse 
                     % Brooks-Corey capillary pressure curve.   Since the
                     % data type is different in the two cases, we cannot set 
                     % the default directly here, but have to set it after
                     % the call to merge_options in the main function.
   
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
                                  reservoirT, pnum, tnum)
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
                                        'reservoirT', reservoirT);
   end
end

% ----------------------------------------------------------------------------

function fun = as_function_of_p(val)

   assert(isnumeric(val) && isscalar(val));
   fun = @(p, varargin) p * 0 + val;

end

% ----------------------------------------------------------------------------
function res = ifelse(cond, yes, no)
    if cond
        res = yes;
    else
        res = no;
    end
end
