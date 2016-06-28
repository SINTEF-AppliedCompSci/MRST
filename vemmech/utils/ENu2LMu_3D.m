function [lambda, mu] = ENu2LMu_3D(E, nu)

   lambda = (E .* nu) ./ ( (1 + nu) .* (1 - 2 * nu));
   
   mu = E ./ (2 * (1 + nu));

end
