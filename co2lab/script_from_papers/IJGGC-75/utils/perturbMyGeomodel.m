function [Gt, rock] = perturbMyGeomodel(Gt, rock, varargin)

    % options:
    opt.perturbType = 'perturbCaprock';
    opt.pert_interval = 5;
    opt = merge_options(opt, varargin{:});
    
    Gt_base = Gt;
    rock_base = rock;
    pert_interval = opt.pert_interval;
    
    
    if strcmpi(opt.perturbType,'perturbCaprock')
        pert_clevel = 1;
        kernal_size = [3 3 3];
        std = 0.65;

    elseif strcmpi(opt.perturbType,'perturbPoro') || strcmpi(opt.perturbType,'perturbPerm')
        kernal_size = [21 3 3];
        std = 3;
        poroRange = opt.pert_interval;
    end
    
    
    if strcmpi(opt.perturbType,'perturbCaprock')
        [Gt,~] = perturb_topSurface(Gt_base, ...
            'pert_interval', pert_interval, ...
            'pert_clevel', pert_clevel, ...
            'kernal_size', kernal_size, ...
            'std', std);
        rock = rock_base;
    
    elseif strcmpi(opt.perturbType,'perturbPoro')
        Gt = Gt_base;
        [rock.poro, ~] = perturb_rockProps(Gt, rock_base.poro, poroRange, ...
            'kernal_size',kernal_size, 'std',std);
        rock.perm = rock_base.perm;

    elseif strcmpi(opt.perturbType,'perturbPerm')
        Gt = Gt_base;
        [rock.poro, ~] = perturb_rockProps(Gt, rock_base.poro, poroRange, ...
            'kernal_size',kernal_size, 'std',std);
        % Function handles for phi-depth, perm-phi model (where
        % coefficients were computed previously using Sto geomodel which
        % contained heterogeneous rock data)
        P_perm = [34.971830212726744 -7.245959920003443];
        permfun = @(poro) exp(P_perm(1) * poro + P_perm(2)); % darcy
        rock.perm = permfun(rock.poro); % darcy
        rock.perm = convertFrom(rock.perm,darcy); % meter^2
        rock.poro = rock_base.poro;
    end

end