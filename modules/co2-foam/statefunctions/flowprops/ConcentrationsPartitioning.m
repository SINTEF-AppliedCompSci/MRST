classdef ConcentrationsPartitioning < StateFunction
    
    methods
        
        %-----------------------------------------------------------------%
        function gp = ConcentrationsPartitioning(model,varargin)
            
            % Parent constructor
            gp@StateFunction(model,varargin{:});
            % Define dependencies
            gp = gp.dependsOn({'Density', 'PoreVolume'}, ...
                'PVTPropertyFunctions');
            gp = gp.dependsOn({'foam', 'water'}, 'state');
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn({'pressure'}, 'state');
            end
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function cm_partitions = evaluateOnDomain(prop, model, state)
            
            fluid = model.fluid;
            G     = model.G;
            
            [cs, sW] = model.getProps(state, 'foam', 'water');
            [rho,pv] = model.getProps(state, 'Density','PoreVolume');
            [rhoW, rhoG] = deal(rho{:});
            
           % Concentration partioning
            if isfield(fluid,'Cpart')
            % Check if field is nonempty
              if ~isempty(fluid.Cpart)
                 % Partitioning constant is given
                 Cpart = fluid.Cpart;
                 assert(Cpart>= 0, 'Partitioning constant must be non-negative.')
              else
                 Cpart = NaN;
              end
            else
              % No partitioning field. foam stays in the phase given by the
              % fluid model, i.e. by fluid.surfingas = true/false.
              Cpart = NaN; % Not used in this case.
            end
	        % Partitioning works as: cG/cW=Cpart, where cG and cW are mass
            % concentrations of surfactant in gas and water, respectively.
	        % Partitioning also needs to be consistent with the adsorption, which
	        % is calculated from the concentration in the water phase, except for
	        % the case where foam is not soluble in water. (This happens if
	        % Cpart is not defined, and fluid.surfingas = true.)
    
            % Resolve various possible configurations
	        [doPart, doAds] = deal(false);
            if ~isnan(Cpart), doPart = true; end
            if isfield(fluid,'adsMax')
               if fluid.adsMax > 0
                  doAds = true;
               end
            end
            mWat = pv.*rhoW.*sW;
            mGas = pv.*rhoG.*(1-sW);
            mRck = (G.cells.volumes-pv).*fluid.rhoRSft;
            % Total mass of foam in each grid cell. The concentration cs is
            % in kg per m3 pore volume. Concentrations cW, cG and cads are
            % in kg per kg.
            mSrf = pv.*cs;
            
           % Consistent set of concentrations conserving foam mass
           [cW, cG, cads] = solvePartitioning(model, fluid, mSrf, ...
               Cpart, mGas, mWat, mRck, doPart, doAds, fluid.surfingas);
       
           cm_partitions{1} = cW;
           cm_partitions{2} = cG;
           cm_partitions{3} = cads;
           cm_partitions{4} = mWat;
           cm_partitions{5} = mGas;
           cm_partitions{6} = mSrf;
           
        end
        %-----------------------------------------------------------------%
        
    end
end

%-------------------------------------------------------------------------%
function [cW, cG, cads] = solvePartitioning(model, fluid, mS, Cpart, mG, mW, mR, doPart, doAds, sInGas)
%Solve partitioning equation for surfactant concentration in water
%
% Input:
% fluid - Fluid model, with adsorption function
% mS    - Total mass of surfactant in each grid cell
% Cpart - Partitioning constant
% mG    - Mass of gas in each grid cell
% mW    - Mass of water in each grid cell
% mR    - Mass of rock in each grid cell (for adsorption)
% doPart
% doAds
% sInGas
%
% Solve the nonlinear equation 
%    mS = mW*cW/(1-cW)+ mG*cG/(1-cG) + ads(cW)*mR 
% in each cell by Newton iterations. cW and cG are the mass 
% concentrations of surfactant in the wate rand gas phases, 
% respectively, and related through the partitioning constant Cpart.
% Variations to the equation in case of no adsorption, no gas
% solubiliy, or no water solubility, are coded individually.

% If mS is an ADI variable, the resulting concentrations get the
% appropriate jacobian matrices added at the end of the function.
    
    if ~doAds
        % No adsorption
        if doPart
            % Partitioning but no adsorption
            % Solve iteratively for cW
            xS = mS;
            mS = value(mS); mW = value(mW); mG = value(mG);
            % Initial approximation when disregarding mass of dissolved 
            % surfactant
            c = mS./(mW+Cpart*mG);
            [c, mS] = model.AutoDiffBackend.initVariablesAD(c, mS);
            func = @(x) mS - mW.*x./(1-x) + Cpart*mG.*x./(1-Cpart*x);
            it = 0;
            maxit = 10;
            abstol = 1.0e-10;
            eq = func(c);
            resnorm = norm(value(eq),'inf');
            while (resnorm > abstol) && (it <= maxit)
                dc = -eq.jac{1}\eq.val;
                c = c + dc;
                c.val(c.val<0) = 0;
                eq = func(c);
                resnorm = norm(value(eq),'inf');
                it = it + 1;
            end
            % Add jacobian information if necessary. Disregard assumed weaker
            % dependency on saturation (water and gas masses).
            if isa(xS, 'ADI')
                dcfdmS = -eq.jac{1}\eq.jac{2};
                cfAD = xS;
                cfAD.val = value(c);
                for i = 1:numel(xS.jac)
                    cfAD.jac{i} = dcfdmS*xS.jac{i};
                end
                cf = cfAD;
            else
                cf = value(c);
            end
            cW = cf;
            cG = Cpart*cW;
        elseif sInGas % No partitioning, surfactant only in gas
            % Safeguard against zero gas saturation. mS should be
            % zero for these cells. 
            % TODO: Consider enforcing this?
            mG = value(mG); mG(mG==0)=1;
            cG = mS./(mS+mG); 
            cW = cG.*0.;
        else % No partitioning, surfactant only in water
            % Safeguard against zero water saturation. mS should be
            % zero for these cells. 
            % TODO: Consider enforcing this?
            mW = value(mW); mW(mW==0)=1;
            cW = mS./(mS+mW);
            cG = cW*0.;
        end
        cads = zeros(size(value(mS)));
    else
        % Have adsorption. Need to find consistent surfactant concentration
        % in the fluid phases.
        % Below c is the concentration in water, except if doPart==false
        % and sInGas==true, in which case c is the concentration in gas.
        xS = mS;
        mS = value(mS); mW = value(mW); mG = value(mG); mR = value(mR);
        scale = max(1./(mS + mW + mG + mR));
        % Start value set to 0 if a significant part of surfactant will be
        % adsorbed, i.e. if mS<mR*fluid.adsMax. Otherwise set start value
        % to mS./mfluids (inside doPArt logic)
        c = mS.*0;
        idx = (mS./mR)>fluid.adsMax;
        if doPart
            c(idx)=mS(idx)./(mW(idx)+Cpart.*mG(idx));
            [c, mS] = model.AutoDiffBackend.initVariablesAD(c, mS);
            func    = @(x) mS - mW.*x./(1-x) + Cpart.*mG.*x./(1-Cpart.*x) - mR.*fluid.surfads(x);
        elseif sInGas % No partitioning, surfactant in gas and adsorbed
            c(idx)=mS(idx)./mG(idx);
            [c, mS] = model.AutoDiffBackend.initVariablesAD(c, mS);
            func    = @(x) mS - mG.*x./(1-x) - mR.*fluid.surfads(x);
        else % No partitioning, surfactant in water and adsorbed
            c(idx)=mS(idx)./mW(idx);
            [c, mS] = model.AutoDiffBackend.initVariablesAD(c, mS);
            func    = @(x) mS - mW.*x./(1-x) - mR.*fluid.surfads(x);
        end
        it      = 0;
        maxit   = 10;
        abstol  = 1e-7;
        eq      = func(c);
        resnorm = norm(value(eq)*scale,'inf');
        is_diag = isa(model.AutoDiffBackend, 'DiagonalAutoDiffBackend');
        while (resnorm > abstol) && (it <= maxit)
            if is_diag
                dc = -eq.val./eq.jac{1}.diagonal(:,1);
            else
                dc = -eq.jac{1}\eq.val;
            end
            c  = c + dc;
            % Restrict concentration to positive values
            c.val(c.val<0) = 0;
            eq = func(c);
            resnorm = norm(value(eq)*scale, 'inf');
            it = it + 1;
        end
        if resnorm > abstol
            warning(['Partitioning calculation did not converge in ', ...
                    '%i iterations. Scaled residual is %.2e'], it - 1, resnorm);
        end
        % Add jacobian information if necessary. Disregard assumed weaker
        % dependency on saturation (water and gas masses).
        if isa(xS, 'ADI')
            if is_diag
                dcfdmS = -eq.jac{1}.diagonal(:,2)./eq.jac{1}.diagonal(:,1);
            else
                dcfdmS = -eq.jac{1}\eq.jac{2};
            end
            cfAD = xS;
            cfAD.val = value(c);
            for i = 1:numel(xS.jac)
                if is_diag
                    cfAD.jac{i}.diagonal = dcfdmS.*xS.jac{i}.diagonal;
                else
                    cfAD.jac{i} = dcfdmS*xS.jac{i};
                end
            end
            cf = cfAD;
        else
            cf = value(c);
        end
        
        if doPart
            cW = cf;
            cG = Cpart*cW;
        elseif sInGas
            cG = cf;
            cW = 0*cG;
        else
            cW = cf;
            cG = 0*cW;
        end
        % Calculate adsorbed surfactant based on concentration in main 
        % fluid phase
        cads = fluid.surfads(cf);
        
    end
    
end
%-------------------------------------------------------------------------%