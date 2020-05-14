classdef WellConMobilityWithPolymer < StateFunction
       
    properties
        % TODO: this is a duplicated property with WellPhaseFlux
        allowCrossFlow = true;
    end
    
    methods
        function gp = WellConMobilityWithPolymer(varargin)
            gp@StateFunction(varargin{:});
            gp = gp.dependsOn({'FacilityWellMapping', 'PerforationPressureGradient', 'WellIndex'});
            gp = gp.dependsOn({'Mobility'}, 'FlowPropertyFunctions');
            gp.label = '\lambda_{wc}'; 
        end
        
        function mobw = evaluateOnDomain(prop, model, state)
            map = prop.getEvaluatedDependencies(state, 'FacilityWellMapping');
            mob = model.ReservoirModel.getProps(state, 'Mobility');
            mobw = cellfun(@(x) x(map.cells), mob, 'UniformOutput', false); 
            
            [dp, wi] = prop.getEvaluatedDependencies(state, 'PerforationPressureGradient', 'WellIndex');
            
            Tdp = -wi.*dp;
            
            mobw = updateMobility(mobw, Tdp, model, state, map, prop); 
        end 
    end   
end


function mobw = updateMobility(mobw, Tdp, model, state, map, prop)
% Currently, this function only handles the effects related to polymer
    vTdp = value(Tdp);
    injection = vTdp > 0;

    % for polymer injecting, we need to modify the injecting
    % mobility
    check = @(prop) isprop(model.ReservoirModel, prop) && model.ReservoirModel.(prop);
    haspolymer = check('polymer');
    if (haspolymer)
        cp = model.ReservoirModel.getProps(state, 'polymer');
        % TODO: we might need to only do this to injection well
        % instead of all the injection perforations
        % TODO: how to handle the crossflow of polymer remains to be
        % fixed
        % TODO: and also, we can trying to have the fully
        % mixing multiplier (or some other thing related) as
        % state function, then a few places related to
        % different viscosities can be built upon these state
        % functions
        wc_inj = (map.cells(injection));
        cp_inj = cp(wc_inj);
        % TODO: this should be the phaseIndex of PolymerComment
        wIx = 1;
        mobw_injw = mobw{wIx}(injection);

        viscpmult = model.ReservoirModel.getProps(state, 'PolymerEffViscMult');
        viscpmult_inj = viscpmult(wc_inj);
        viscpmultfull = model.ReservoirModel.fluid.muWMult(cp_inj);

        mobw{wIx}(injection) = mobw_injw ./ viscpmultfull .* viscpmult_inj;
                
        % the shear effects for polymer
        hasShear= check('usingShear') || check('usingShearLog') || check('usingShearLogshrate');
        if hasShear
            q_ph = calculatePhaseRate(Tdp, mobw, map, prop.allowCrossFlow, model);
            mobw = applyShearEffectsWell(mobw, q_ph, prop, model.ReservoirModel, state);
        end
    end
end