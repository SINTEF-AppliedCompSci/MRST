function [Ieta, varargout] = etaIntegrals(EOS, P, T, Gct, gct)
% 
% SYNOPSIS:
%   function [Ieta, varargout] = etaIntegrals(EOS, P, T, Gct, gct)
%
% PARAMETERS:
%   EOS - Equation of state, containing the following density-related
%         functions of P and T: 'rho', 'beta', 'gamma', 'beta2', 'gamma2',
%         'chi'
%   P   - reference pressure                         (scalar or vector)
%   T   - reference temperature                      (scalar or vector)
%   Gct - vertical component of temperature gradient (scalar, deg/m)
%   gct - vertical component of gravity vector       (scalar)
%
% RETURNS:
%   Ieta      - Function of one argument 'h', giving the approximated
%               integral of the 'eta'-function  from reference to reference+h.
%   varargout - INupEta - Integral of 'eta' multiplied by Nu_p.
%               INugEta - Integral of 'eta' multiplied by Nu_G.
%               Ieta2   - Integral of 'eta' squared.
%               Eta     - (non -integrated) eta -function
%               Nup     - (non -integrated) Nu_p function (dP/dP_ref)    
%               Nug     - (non -integrated) Nu_G function (dP/dT)
% EXAMPLE:
%
% SEE ALSO: eta, 
%
    is_full = (strcmpi(EOS.compressible, 'full'));
    
    if ~is_full
        if (isa(P,'ADI'))
            dim = size(P.val);
        else
            dim = size(P);
        end
        rho = EOS.rho(P, T);
        Ieta = @(h) ones(dim);
        if (nargout > 1) varargout{1} = @(h) ones(dim)  ; end ;  % INupEta
        if (nargout > 2) varargout{2} = @(h) zeros(dim) ; end ;  % INugEta
        if (nargout > 3) varargout{3} = @(h) ones(dim)  ; end ;  % Ieta2
        if (nargout > 4) varargout{4} = @(h) ones(dim)  ; end ;  % Eta
        if (nargout > 5) varargout{5} = @(h) ones(dim)  ; end ;  % Nup
        if (nargout > 6) varargout{6} = @(h) zeros(dim) ; end ;  % Nug
        return;
    end

    % If we got here, we are in the fully compressible case
    rho_g  = EOS.rho    (P, T) .* gct;
    beta   = EOS.beta   (P, T);
    gamma  = EOS.gamma  (P, T);
    beta2  = EOS.beta2  (P, T);
    gamma2 = EOS.gamma2 (P, T);
    chi    = EOS.chi    (P, T);

    Ieta = approximateIeta(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    if (nargout > 1) 
        varargout{1} = approximateINupEta(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    end;
    if (nargout > 2)
        varargout{2} = approximateINugEta(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    end;
    if (nargout > 3) 
        varargout{3} = approximateIEta2 (rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    end;
    if (nargout > 4)
        varargout{4} = eta(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    end
    if (nargout > 5)
        varargout{5} = Nup(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    end
    if (nargout > 6)
        varargout{6} = Nug(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    end
end

function fct = Nup(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
    [d dd] = NupDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    fct = taylorDev(1, d, dd);
end

function fct = Nug(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
    [d dd] = NugDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    fct = taylorDev(0, d, dd);
end


function fct = eta(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
    [d dd] = etaDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    fct = taylorDev(1, d, dd);
end

function fct = taylorDev(f0, fd, fdd)
    fct = @(h) f0 + (fd.*h) + 0.5 * (fdd .* h .* h);
end

function [d dd] = etaDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
    
    d  = (rho_g .* beta) - gamma .* Gct;
    dd = ((rho_g.^2) .* (beta.^2 + beta2)) + ...
         (rho_g .* Gct .* ((2 .* chi) - (gamma .* beta))) + ...
         Gct.^2 .* gamma2;
end

function [d dd] = NupDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct)

    d = rho_g .* beta;
    dd = rho_g .* ((rho_g .* (beta.^2 + beta2)) + (Gct .* chi));

end

function [d dd] = NugDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
    
    d =  -rho_g .* gamma;
    dd =  rho_g .* (rho_g .* (chi - gamma.*beta) + Gct.*gamma2);
end


function Ieta = approximateIeta(rho_g, beta, gamma, beta2, gamma2, chi, Gct)

    [d, dd] = etaDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    Ieta = taylorDev(1, 1/2 * d, 1/3 * dd);
end

function INupEta = approximateINupEta(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
        
    [eta_d, eta_dd] = etaDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    [nup_d,  nup_dd ] = NupDer (rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    
    d = eta_d + nup_d;
    dd = eta_dd + 2 .* eta_d .* nup_d + nup_dd;
    
    INupEta = taylorDev(1, 1/2 * d, 1/3 * dd);
end

function INugEta = approximateINugEta(rho_g, beta, gamma, beta2, gamma2, chi, Gct)

    [eta_d, eta_dd] = etaDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    [nug_d, nug_dd] = NugDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);

    d = nug_d;  % eta' * nu + eta * nu' , with eta == 1 and nu == 0
    dd = nug_dd + 2 * eta_d .* nug_d;
    
    INugEta = taylorDev(0, 1/2 * d, 1/3 * dd);
    
end

function IEta2 = approximateIEta2(rho_g, beta, gamma, beta2, gamma2, chi, Gct)
    
    [eta_d, eta_dd] = etaDer(rho_g, beta, gamma, beta2, gamma2, chi, Gct);
    
    d  = 2 * eta_d;
    dd = 2 * (eta_dd + eta_d.^2);
    
    IEta2 = taylorDev(1, 1/2 * d, 1/3 * dd);
end

