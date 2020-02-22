classdef PhasePotentialDifferenceMPFA < StateFunction
    % The difference in potential over a face as defined for a MPFA scheme
    properties
        ViscousGrad % Matrix mapping cell pressures to interior face phase gradient
        GravityGrad % Matrix mapping the interior faceheight differences to potential differences
        uniquePressure; % Is there a single unique pressure? If so, we can speed up the calculations
        hasGravity; % Is there gravity?
    end
    
    methods
        function pp = PhasePotentialDifferenceMPFA(model, varargin)
            pp@StateFunction(model);
            if isempty(model.FlowPropertyFunctions)
                hasPC = true; % We do not know. Take the performance hit.
            else
                hasPC = model.FlowPropertyFunctions.CapillaryPressure.pcPresent(model);
            end
            pp.uniquePressure = ~hasPC;
            pp.hasGravity = norm(model.gravity(1:model.G.griddim)) > 0;
            require mpfa
            [~, M] = computeMultiPointTrans(model.G, model.rock);
            
            Tv = M.rTrans;
            Tg = M.rgTrans(model.operators.internalConn, :);
            
            scale = -1./(2*model.operators.T);
            pp.ViscousGrad = Tv.*scale;
            pp.GravityGrad = -Tg.*scale/2;
            assert(all(M.N(:) == model.operators.N(:)));
            
            pp = merge_options(pp, varargin{:});
            pp.label = 'K(\nabla p_\alpha+g\rho_\alpha\Delta z)';
        end
        
        function v = evaluateOnDomain(prop, model, state)
            nph = model.getNumberOfPhases();
            Vgrad = prop.ViscousGrad;
            Ggrad = prop.GravityGrad;
            v = cell(1, nph);
            if prop.uniquePressure
                p = model.getProp(state, 'pressure');
                Kgradp = Vgrad*p;
                [v{:}] = deal(Kgradp);
            else
                pressures = model.getProp(state, 'PhasePressures');
                for i = 1:nph
                    v{i} = Vgrad*pressures{i};
                end
            end
            
            if prop.hasGravity
                gdz = model.getProp(state, 'GravityPotentialDifference');
                for i = 1:nph
                    v{i} = v{i} + Ggrad*gdz{i};
                    % v{i} = v{i} + gdz{i}; Alternate way, inconsistent
                end
            end
        end
    end
end