function fluid = makeADIFluid(G, type, varargin)
%
% Construct an ADI fluid with properties specific to a chosen model
% 
% SYNOPSIS:
%   function fluid = makeADIFluid(type, opt, varargin)
%
% DESCRIPTION:
%
% @@ _complete_me_
%
% PARAMETERS:
%   G            - @@ 
%   type         - @@
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

   %mu(3), rho(3), n(3), sr, sw, pvMultR, bW, bG, surface_tension
   opt.rhoS = ; % reference brine and gas densities
   opt.mu   = ; % brine and gas viscosities (constants or functions)
   opt.res  = ; % residual brine and residual gas

   %% adding type-specific modifications
   switch type
     case 'simple' %@@ tested anywhere?
       fluid = makeSimpleFluid(G, opt);
     case 'integrated' %@@ tested anywhere? 
       fluid = makeIntegratedFluid(G, opt);
     case 'sharp interface'
       fluid = makeSharpInterfaceFluid(G, opt);
     case 'linear cap.'
       fluid = makeLinCapFluid(G, opt);     
     case 'S table'
       fluid = makeSTableFluid(G, opt);
     case 'P-scaled table'
       fluid = makePScaledFluid(G, opt);
     case 'P-K-scaled table'
       fluid = makePKScaledFluid(G, opt);
     otherwise
       error([type, ': no such fluid case.']);
   end
         
end

% ============================================================================

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

function fluid = constant_density(fluid, rhoS)
   
   fluid = setfield(fluid, 'rhoWS', rhoS(1));
   fluid = setfield(fluid, 'rhoGS', rhoS(2));
   
end

% ----------------------------------------------------------------------------

function fluid = constant_viscosity(fluid, mu)

   fluid = setfield(fluid, 'muW', as_function_of_p(mu(1)));
   fluid = setfield(fluid, 'muG', as_function_of_p(mu(2)));
   
end

% ----------------------------------------------------------------------------

function fluid = linear_relperms(fluid)

   fluid = setfield(fluid, 'krW' , @(sw, varargin) sg);
   fluid = setfield(fluid, 'krG' , @(sg, varargin) sg);
   fluid = setfield(fluid, 'kr3D', @(s           )  s);
   
end

% ----------------------------------------------------------------------------

function fluid = residual_saturations(fluid, resSat)

   fluid = setfield(fluid, 'res_water', resSat(1));
   fluid = setfield(fluid, 'res_gas'  , resSat(2));
   
end

% ----------------------------------------------------------------------------

function fluid = constant_formation_factors(fluid)

   fluid = setfield(fluid, 'bW', as_function_of_p(1));
   fluid = setfield(fluid, 'BW', as_function_of_p(1));
   fluid = setfield(fluid, 'bG', as_function_of_p(1));
   fluid = setfield(fluid, 'BG', as_function_of_p(1));
   
end

% ----------------------------------------------------------------------------

function fluid = sharp_interface_cap_pressure(fluid, G)
   
   % When function is called, the following fields are needed
   require_fields(fluid, {'bW', 'bG', 'rhoWS', 'rhoGS'});
   
   fluid = setfield(fluid, 'pcWG', @(sg, p, varargin)                       ...
                                    norm(gravity) *                         ...
                                    (fluid.rhoWS .* fluid.bW(p) -           ...
                                     fluid.rhoGS .* fluid.bG(p)) .* (sg) .* ...
                                    G.cells.H); 
   
   fluid = setfield(fluid, 'invPc3D', @(p) 1 - (sign(p + eps) + 1) / 2);

end

% ============================================================================

function fluid = makeSimpleFluid(G, opt)

   fluid = constant_density([], opt.rhoS);          % 'rhoWS'    , 'rhoGS'
   fluid = constant_viscosity(fluid, opt.mu);       % 'muW'      , 'muG'
   fluid = linear_relperms(fluid);                  % 'krW'      , 'krG', 'kr3D'
   fluid = residual_saturations(fluid, opt.resSat); % 'res_water', 'res_gas'
   fluid = constant_formation_factors(fluid);       % 'bW', 'BW' ,  'bG', 'BG'
   fluid = sharp_interface_cap_pressure(fluid, G);  % 'pcWG'     , 'invPc3D'
   
end

% ----------------------------------------------------------------------------

function fluid = makeIntegratedFluid(G, opt)
   
   fluid = addVERelpermIntegratedFluid;
   

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
   
   
   
end

% ----------------------------------------------------------------------------

function fluid = makeSharpInterfaceFluid(G, opt)
   
end

% ----------------------------------------------------------------------------

function fluid = makeLinCapFluid(G, opt)
   
end

% ----------------------------------------------------------------------------

function fluid = makeSTableFluid(G, opt)
   
end

% ----------------------------------------------------------------------------

function fluid = makePScaledFluid(G, opt)
   
end

% ----------------------------------------------------------------------------

function fluid = makePKScaledFluid(G, opt)
   
end

% ----------------------------------------------------------------------------

