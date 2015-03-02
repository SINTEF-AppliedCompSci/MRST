function obj = CO2props_new(varargin)
   
   opt.rhofile              = 'rho_huge';
   opt.mufile               = 'mu_huge';
   opt.hfile                = 'h_small';
   opt.const_derivatives    = false;
   opt.assert               = true;
   opt.sharp_phase_boundary = true;
   opt = merge_options(opt, varargin{:});
   
   boundary = [];
   if opt.sharp_phase_boundary
      [p_c, t_c]= CO2CriticalPoint();
      boundary = {[p_c, t_c], @CO2VaporPressure};
   end

   setupfun = @(name, file) ...
       SampledProp2D(name, file, ...
                     'phase_boundary'   , boundary             , ...
                     'const_derivatives', opt.const_derivatives, ...
                     'assert_in_range'  , opt.assert);

   obj = setupfun('rho', opt.rhofile);
   obj = add_fields(obj, setupfun('mu' , opt.mufile));
   obj = add_fields(obj, setupfun('h'  , opt.hfile));

      
   %%%%
   obj.beta  = @(P, T)   obj.rhoDP(P, T) ./ obj.rho(P, T);  % Compressibility coef.
   obj.gamma = @(P, T) - obj.rhoDT(P, T) ./ obj.rho(P, T);  % Coef. of thermal expansion
  
  % compressibility factor cross derivative (equal for both Pbeta and Tbeta)
  obj.betaDPDT = @(P,T) (obj.rhoDPT(P, T) + obj.rhoDP(P, T) .* obj.rhoDT(P, T)) ./ obj.rho(P, T);

end

% ----------------------------------------------------------------------------

function obj = add_fields(obj, other)

   to_add = fields(other);
   for i = 1:numel(to_add)
      obj.(to_add{i}) = other.(to_add{i});
   end
end
