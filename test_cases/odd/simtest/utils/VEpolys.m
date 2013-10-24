function polys = VEpolys

polys.rgh        = @rgh; % pstar
polys.alpha      = @alpha_poly;
polys.intAlpha   = @intAlpha_poly; %hpi
polys.intSqAlpha = @intSqAlpha_poly; % A

end

%% Implementation of member functions

function p = rgh(h, rho_top, theta)
  p = norm(gravity) * cos(theta) * rho_top .* h;
end

function res = alpha_poly(h, rho_top, beta_top, beta_der_top, theta)
  p = rgh(h, rho_top, theta);
  res = 1 + ...
        beta_top .* p + ...
        1/2 * (2 * beta_top.^2 + beta_der_top) .* p.^2;
end

function res = intAlpha_poly(h, rho_top, beta_top, beta_der_top, theta)
  p = rgh(h, rho_top, theta);
  res = 1 + ...
        1/2 * beta_top .* p + ...
        1/6 * (2 * beta_top.^2 + beta_der_top) .* p.^2;

  res = res .* h;
end


function res = intSqAlpha_poly(h, rho_top, beta_top, beta_der_top, theta)
% expresses the definite integral of 'alpha' squared, from top to 'h'.
  p = rgh(h, rho_top, theta);
  
  % we truncate terms higher than p.^2
  res = 1 + ...
        beta_top .* p + ...
        (beta_top.^2 + (1/3 * beta_der_top)) .* p.^2; 
  
  res = res .* h;
  
end

