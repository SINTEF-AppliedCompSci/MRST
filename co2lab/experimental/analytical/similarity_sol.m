function [pfun, hfun, zones] = similarity_sol(H, phi, k, Qvol, sb_res, k_rc, ...
                                              mu_c, mu_b, ...
                                              rho_c, rho_b, c_r, c_b)
%% Intermediate definitions
    lambda_c = k_rc / mu_c;
    lambda_b = 1/mu_b;
    drho = rho_b - rho_c;
    g = 9.81; % Gravitational constant

%% Defining dimensionless groups
    lambda = lambda_c / lambda_b;
    Gamma  = (2*pi*drho*g*k*lambda_b*H^2)  / (Qvol);
    Psi    = (35.33*pi*phi*H*k*lambda_b*(1-sb_res)) / ((c_r + phi*c_b)*Qvol);
        
%% Defining helper functions
    chi_const_fac = 2*pi*phi*H*(1-sb_res)/Qvol;
    chi_fun   = @(r, t) chi_const_fac * r.^2 ./ t;
    hprim_fun = @(chi) h_impl(chi, lambda);
    F_fun     = @(h) -lambda/(lambda-1) * (h - log((lambda-1)*h +1) / (lambda-1));
    
    Pdb_fun = @(chi) Pdb_fun_impl(chi, lambda, Psi, Gamma, F_fun(hprim_fun(chi)));
    
%% Result functions
    hfun = @(r, t) H * hprim_fun(chi_fun(r, t));
    pfun = @(r, t, r0, t0, p0) p0 + ...
                               drho*g*H * (Pdb_fun(chi_fun(r, t)) - ...
                                           Pdb_fun(chi_fun(r0, t0)));
%% Result zone markers
    
    zones = @(t) sqrt([2/lambda, 2*lambda, Psi] .* t ./ chi_const_fac);
end

%% HELPER FUNCTIONS

% Implementation of function computing dimensionless height
function h = h_impl(chi, lambda)
    z1_ix = (chi <= 2/lambda);
    z2_ix = (chi <  2*lambda) & ~z1_ix;
    z3_ix = (chi >= 2*lambda);
    
    h = zeros(size(chi));
    h(z1_ix) = 1;
    h(z2_ix) = inv(lambda-1) * (sqrt(2*lambda./chi(z2_ix))-1) ;
    h(z3_ix) = 0;
end

function dp = Pdb_fun_impl(chi, lambda, Psi, Gamma, F)
    assert(Psi > 2*lambda);
    
    dp = zeros(size(chi));

    % Specification of delta-pressure functions
    dp_z4 = @(x) zeros(size(x));
    dp_z3 = @(x) 1/(2*Gamma)*W(35.33*x/(8*Psi)) + ...
                 dp_z4(Psi*ones(size(x))) ;
    dp_z2 = @(x) 1/Gamma - sqrt(x)/(Gamma*sqrt(2*lambda)) + ...
                 dp_z3(2*lambda*ones(size(x)));
    dp_z1 = @(x) -1/(2*lambda*Gamma) * log(x*lambda/2) + ...
                 dp_z2(2/lambda*ones(size(x)));
    
    z4_ix = (chi >= Psi);
    z3_ix = (chi >= 2*lambda) & ~(z4_ix);
    z2_ix = (chi >= 2/lambda) & ~(z4_ix | z3_ix);
    z1_ix = ~(z4_ix | z3_ix | z2_ix);
    
    % Computing delta-pressure for the values of chi
    dp(z4_ix) = dp_z4(chi(z4_ix));
    dp(z3_ix) = dp_z3(chi(z3_ix));
    dp(z2_ix) = dp_z2(chi(z2_ix));
    dp(z1_ix) = dp_z1(chi(z1_ix));
    
    % Correction to get bottom-pressure
    dp = dp + F;
end



% function dp = Pdb_fun_impl(chi, lambda, Psi, Gamma, F_fun, h_fun)
%     assert(Psi > 2*lambda);
    
%     dp = zeros(size(chi));
    
%     dp_z4 = @(x) zeros(size(x)) + F_fun(h_fun(Psi*ones(size(x))));
%     dp_z3 = @(x) 1/(2*Gamma)*W(35.33*x/(8*Psi)) + ...
%                  z4_fun(Psi*ones(size(x))) + ...
%                  F_fun(h_fun(x));
%     dp_z2 = @(x) 1/Gamma - sqrt(x)/(Gamma*sqrt(2*lambda)) + ...
%                  z3_fun(2*lambda*ones(size(x))) + ...
%                  F_fun(h_fun(x));
%     dp_z1 = @(x) -1/(2*lambda*Gamma) * log(x*lambda/2) + ...
%                  z2_fun(2/lambda*ones(size(x))) + ...
%                  F_fun(h_fun(x));
    
%     z4_ix = (chi >= Psi);
%     z3_ix = (chi >= 2*lambda) & ~(z4_ix);
%     z2_ix = (chi >= 2/lambda) & ~(z4_ix | z3_ix);
%     z1_ix = ~(z4_ix | z3_ix | z2_ix);
    
%     dp(z4_ix) = dp_z4(chi(z4_ix));
%     dp(z3_ix) = dp_z3(chi(z3_ix));
%     dp(z2_ix) = dp_z3(chi(z2_ix));
%     dp(z1_ix) = dp_z1(chi(z1_ix));
% end


    
    