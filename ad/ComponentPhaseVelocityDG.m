classdef ComponentPhaseVelocityDG < StateFunction
    properties (Access = protected)
        mobility_name; % Name of state function where component mobility comes from
    end

    methods
        function cf = ComponentPhaseVelocityDG(backend, mob_name)
            cf@StateFunction(backend);
            if nargin < 2
                mob_name = 'ComponentMobility'; 
            end
            cf.mobility_name = mob_name;
            cf = cf.dependsOn({'PressureGradient'});
            cf = cf.dependsOn({cf.mobility_name, 'GravityPermeabilityGradient'}, 'FlowPropertyFunctions');
        end

        function v = evaluateOnDomain(prop, model, state)
            ncomp   = model.getNumberOfComponents;
            nph     = model.getNumberOfPhases;
            dp      = prop.getEvaluatedDependencies(state, 'PressureGradient');
            gRhoKdz = model.getProps(state, 'GravityPermeabilityGradient');
            compMob = model.getProps(state, prop.mobility_name);
            
            dpc = cell(1, nph);
            ic = model.operators.internalConn;
            [~, n] = rlencode(state.faces);
            n = rldecode(n,n,1);
            i2a = sparse(state.faces, 1:numel(state.faces), 1./n, model.G.faces.num, numel(state.faces));
            for i = 1:nph
                dpc{i} = i2a*dp{i};
                dpc{i} = dpc{i}(state.cells);
            end
            
            v = cell(ncomp, nph);
            for c = 1:ncomp
                for ph = 1:nph
                    mob = compMob{c, ph};
                    if ~isempty(mob)
                        v{c, ph} = -mob.*(dpc{ph} - gRhoKdz{ph});
                    end
                end
            end
        end
    end
end