function bG = boCO2(T_ref, rhoG, varargin)
    
    opt.rho_datafile = 'rho_big_trunc';
    opt.sharp_phase_boundary = true;
    opt = merge_options(opt, varargin{:});
    
    obj= CO2props('rhofile', opt.rho_datafile, ...
                  'sharp_phase_boundary', opt.sharp_phase_boundary);
    rhoCO2 =@(p) obj.rho(p, T_ref);
    bG =@(p) rhoCO2(p)/rhoG;
end