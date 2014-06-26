function bG  =boCO2(T_ref, rhoG, varargin)
    
    opt.datafile = 'rho_big_trunc';
    opt.sharp_phase_boundary = true;
    opt = merge_options(opt, varargin{:});
    
% @@ Warning: this function will cause a memory leak, since the `dispose`
% function is never called on the created `CO2props` object.  (A copy of the
% tables used for interpolation will remain indefinitely in memory).
      %obj= CO2props('rho_small', 'h_small');
      obj= CO2props('rhofile', opt.datafile, ...
                    'sharp_phase_boundary', opt.sharp_phase_boundary);
      rhoCO2 =@(p) obj.rho(p, T_ref);
      %rhoG = rhoCO2(p_ref);
      bG =@(p) rhoCO2(p)/rhoG;
end