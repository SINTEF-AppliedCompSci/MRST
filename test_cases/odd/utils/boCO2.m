function [bG rhoG] =boCO2(p_ref, T_ref)
% @@ Warning: this function will cause a memory leak, since the `dispose`
% function is never called on the created `CO2props` object.  (A copy of the
% tables used for interpolation will remain indefinitely in memory).
      obj= CO2props();
      rhoCO2 =@(p) obj.rho(p, T_ref);
      rhoG = rhoCO2(p_ref);
      bG =@(p) rhoCO2(p)/rhoG;
end