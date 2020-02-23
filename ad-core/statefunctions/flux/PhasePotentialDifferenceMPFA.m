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
            Tv = M.rTrans; % Cells -> Inner faces
            Tg = M.rgTrans(model.operators.internalConn, :); % Inner -> Inner
            % Change sign and re-scale operators to fit with AD-OO
            % definition of gradient.
            scale = -1./(2*model.operators.T);
            pp.ViscousGrad = Tv.*scale;
            pp.GravityGrad = -Tg.*scale/2;
            assert(all(M.N(:) == model.operators.N(:)));
            % Allow kw inputs
            pp = merge_options(pp, varargin{:});
            pp.label = '\Theta_\alpha';
        end
        
        function v = evaluateOnDomain(prop, model, state)
            nph = model.getNumberOfPhases();
            M_d = prop.ViscousGrad;
            M_g = prop.GravityGrad;
            v = cell(1, nph);
            if prop.uniquePressure
                p = model.getProp(state, 'pressure');
                gradp = M_d*p;
                [v{:}] = deal(gradp);
            else
                pressures = model.getProp(state, 'PhasePressures');
                for i = 1:nph
                    v{i} = M_d*pressures{i};
                end
            end
            
            if prop.hasGravity
                gdz = model.getProp(state, 'GravityPotentialDifference');
                for i = 1:nph
                    v{i} = v{i} + M_g*gdz{i};
                    % v{i} = v{i} + gdz{i}; Alternate way, inconsistent
                end
            end
        end
    end
end