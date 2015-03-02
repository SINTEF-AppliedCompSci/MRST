function [P FpEta Int_FpEta stopped odeObj] = numInt_P_FpEta(zref, Pref, Tref, Tgrad, h, fluidProps)
%
% Compute pressure ('P'), \eta x f_p ('FpEta') at ('zref' + 'h'), as well as
% the integrated value of 'FpEta' from 'zref' to 'zref + h'.
%
% \eta(z) is defined as the ratio of density at depth(z) as compared with
% reference density at depth 'zref'.
% f_p  is defined as the partial derivative of pressure at (z) with regards
% to the pressure at the reference depth ('Pref').
% 
% The significance of 'FpEta' is the ratio of the mass flux magnitude at
% depth 'zref + h' over the mass flux at reference depth 'zref', assuming
% that the caprock is level.  Thus, the value 'Int_FpEta' is the ratio of
% total mass flux over the mass flux that would have taken place in an
% incompressible setting (all while assuming a flat caprock).
%
% SYNOPSIS:
%   function [P FpEta Int_FpEta] = numInt_P_FpEta(zref, Pref, h, rho_fun, T_fun) 
%
% PARAMETERS:
%   zref       - reference depth
%   Pref       - pressure at reference depth
%   Tref       - temperature at reference depth
%   Tgrad      - temperature gradient (considered constant), in BASE UNITS (deg/meter)    
%   h          - depth down to which we want to integrate (i.e. to 'zref' + 'h')
%   fluidProps - fluid property objects, providing the following functions
%                * rho(P, T)
%                * rhoDP(P, T)
%                * rhoDT(P, T)
%
% RETURNS:
%   P         - pressure at depth 'zref' + 'h'
%   FpEta     - value of the product of \eta and f_p at depth 'zref' + 'h'
%   Int_FpEta - depth integrated value of 'FpEta' from 'zref' to 'zref' + 'h'.
%   stopped   - true if the integration was prematurely stopped (liqid-vapor boundary hit)
%   odeObj    - the full object returned from the ODE solver (in case it's desired)
%
% EXAMPLE:
%
% SEE ALSO:
%

    %% Defining the system
    g = norm(gravity);
    v0 = [Pref 1 0]';
    
    % The unknown vector v = [pressure, FpEta, int_FpEta]
    function vder = dfun(z, v)
        T     = Tref + (z-zref) * Tgrad;
        P     = v(1);
        rho   = fluidProps.rho(P, T);
        beta  = fluidProps.rhoDP(P, T)   / rho;
        gamma = - fluidProps.rhoDT(P, T) / rho;
        
        vder = [rho*g; ...                                  % dP/dz
                (2*rho*g*beta - gamma * Tgrad) * v(2); ...  % d(FpEta)/dz
                v(2)];                                      % d(int(FpEta, z))/dz = FpEta
    end

    
    % event function that triggers if we cross the liquid-vapor boundary
    last_phase = -1; % 0 = supercritical, 1 = liquid, 2 = gas; -1 = uninitialized
    function [val, isterminal, dir] = trigger(z, v)
        dir = 0;           % We do not care about this variable, but the solver needs it.
        isterminal = true; % Terminate if boundary crossing is detected.
        cur_phase = fluidProps.phaseOf(v(1), Tref + (z-zref) * Tgrad);
        if (last_phase + cur_phase == 3)
            val = 0;  % we have passed from liquid to gas, or vice versa
        else
            val = 1;
        end
        last_phase = cur_phase;
    end

    %% solving the system
    
    %v = ode23(@dfun, [zref, zref+h], v0);
    v = ode23(@dfun, [zref, zref+h], v0, odeset('Events', @trigger));
    
    %% Unpacking result and returning it

    P         = v.y(1, end);
    FpEta     = v.y(2, end);
    Int_FpEta = v.y(3, end);
    if nargout > 3
        stopped = isfield(v, 'xe') && (numel(v.xe) > 0);
        % if stopped
        %     'stopped'
        % end            
    end
    
    if nargout > 4
        odeObj = v;
    end
end

