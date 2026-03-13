classdef MultipliedMobility < Mobility 
% Mobility with foam-dependent multiplier for the gas phase

    methods
        %-----------------------------------------------------------------%
        function pp = MultipliedMobility(model, varargin)
            
            % Call parent constructor
            pp@Mobility(model, varargin{:});
            % Define dependencies
            pp = pp.dependsOn('ConcentrationsPartitioning');
            pp = pp.dependsOn('PermeabilityPotentialGradient', ...
                'FlowDiscretization');
            pp = pp.dependsOn({'foam'}, 'state');
            if isfield(model.fluid, ...
                    {'tranMultR', 'transMult', 'mobMultIsLinear'})
                pp = pp.dependsOn({'pressure'}, 'state');
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function mob = evaluateOnDomain(prop, model, state)
        % Get effective pore-volume, accounting for rock compressibility

            fluid = model.fluid;
            
            mob  = evaluateOnDomain@Mobility(prop, model, state);
            mobG = mob{2};

            % Transmissibility times Phase potentials
            kgrad = model.getProp(state, 'PermeabilityPotentialGradient');
            Cwga  = prop.getEvaluatedDependencies(state,'ConcentrationsPartitioning');
            sW    = model.getProp(state, 'water');
            cW    = Cwga{1};
            
            vGm = []; % Absoulte value of gas velocity at cell centroids
            
            % Map kgrad for gas from faces to cells
            TdpGm = model.operators.sqVeloc(kgrad{2}).^0.5;
            
            sWm = sW;
            % TODO:
            % This should be cG for cases where surfactant is only soluble
            % in gas. I.e. when no partitioning and surfingas==true.
            % For now assume surfactant in water, and use cW.
            cm = cW;

            if fluid.mobMultIsLinear
                [Fcap, Fgrad, Gc] = fluid.mobMult(sWm, cm);
                if fluid.noShear
                    vGm = ones(size(cW));
                else
                    if isempty(vGm)
                        vGm = Fcap.*mobG.*TdpGm./(1 - Fgrad.*mobG.*TdpGm);
                    end
                end
                F = Fcap + Fgrad.*vGm;
                % Add concentration effect
                F = 1 - Gc.*(1-F);
            else
                %if isempty(vGm)
                % Solve nonlinear problem to get vGm
                %
                % Initialize gas velocity using foam strength at reference velocity
                vGm0 = fluid.mobMult(sWm, cm, fluid.foam.v0).*mobG.*TdpGm;
                % vGm0 = face2cell(VG, model);
                vGm  = solveMobMult(fluid, sWm, cm, mobG, TdpGm, vGm0);
                %end
                F = fluid.mobMult(sWm, cm, vGm);
            end
            if fluid.usePermDep
                PermDep = fluid.foam.PermDep(model.rock.perm(:,1)); % Evaluate on x-permeability
                PermDep(PermDep<eps)=eps;
                cmin = 1e-10; cmax = 1e-5;
                f = (value(cm) - cmin)./(cmax - cmin);
                f = max(min(f,1),0);
                PermDep = PermDep.*f + (1-f);
                F = F./PermDep;
            end

            sF = fluid.foam.sF;
            F = min(F,1).*(sW >= sF) + (sW < sF);

            mob{2} = F.*mobG;

        end
        
    end
    %-----------------------------------------------------------------%
    
end

%-------------------------------------------------------------------------%
function vG = solveMobMult(fluid, sW, c, mobG, TdpG, vG0)
% Solve nonlinear equation for gas velocity, when foam mobility effect is 
% velocity dependent.
%
% Parameters:
%   fluid - Fluid definition. Must contain mobiltiy multiplier function fluid.mobMult
%   sW    - Water saturation. Input to mobMult
%   c     - Surfactant concentration. Input to mobMult
%   mobG  - Gas mobility without foam effect. (Black-oil mobility)
%   TdpG  - Gas potential.
%   vG0   - Initial guess for gas velocity.
%
% Returns:
%   vG    - Calculated gas velocity when foam effect on mobility is taken into account.

    % Save (potential) ADI variable with derivatives information.
    x = TdpG;
    [vG, TdpG] = initVariablesADI(value(vG0), value(TdpG));
    
    % Residual equation. We use |vG| and |grad p|, hence the minus sign
    sW = value(sW); c = value(c); mobG = value(mobG);
    
    F = @(v) fluid.mobMult(sW, c, v);
    func = @(v) v - F(v).*mobG.*TdpG;
    
    it = 0;
    maxit = 30;
    abstol = 1.e-10;

    eq = func(vG);
    resnorm = norm(value(eq), 'inf');
    
    while (resnorm > abstol) && (it <= maxit)

        dvG = -eq.jac{1}\eq.val;
        vG = vG + dvG;
        eq = func(vG);
        resnorm = norm(value(eq), 'inf');
    
        it = it + 1;

    end
    
    % Use chain rule to calculate jacobians
    if isa(x, 'ADI')
        dvdTdpG = -eq.jac{1}\eq.jac{2}; 
        vGAD = x;
        vGAD.val = value(vG);
        for i = 1:numel(x.jac)
            vGAD.jac{i} = dvdTdpG*x.jac{i};
        end
        vG = vGAD;
    else
        vG = value(vG);
    end

    if (it >= maxit) && (resnorm > abstol)
        error(['Convergence failure within %d iterations\n', ...
                'Final residual = %.8e'], maxit, resnorm);
    end

end
%-------------------------------------------------------------------------%