function p = analyticmandel(t, x, params, varargin)
% Compute analytical solution for pressure in normalized coordinates
    
    opt = struct('num_modes', 100);
    opt = merge_options(opt, varargin{:});
    
    N = opt.num_modes;
    
    nu = params.nu;
    iM = params.iM;
    
    % only incompressible case for the moment
    assert(iM == 0, 'only consider incompressible case for the moment');
    eta = (1 - nu)/(1 - 2*nu);
    
    % compute values of alphas
    alphas = computealphas(eta, N);

    docheck = false;
    if docheck
        % check alphas
        figure
        hold on
        a = (0 : 1e-2: 10)';
        plot(a/pi, tan(a)/pi);
        plot(a/pi, 2*eta*a/pi);
        axis([0, 10/pi, 0, 2*eta*10/pi]);
    end
    
    p = 0*x;
    casename = 'coussy';
    switch casename
      case 'coussy'
        for i = 1 : N
            % formula from Coussy
            an = alphas(i);
            cn = cos(an*x) - cos(an);
            cn = (sin(an)/(an - sin(an)*cos(an)))*cn;
            cn = exp(-an^2*t)*cn;
            p = p + 2*cn;
        end
      case 'verruijt'
        for i = 1 : N
            % formula from Verruijt (We have checked that they are equivalent to Coussy - up a constant coefficient)
            an = alphas(i);
            cn = cos(an*x) - cos(an);
            cn = cos(an)/(1 - 2*eta*(cos(an).^2))*cn;
            cn = cn*exp(-an^2*t);
            p = p + 2*eta*cn;
        end
    end

end
