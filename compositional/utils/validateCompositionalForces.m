function forces = validateCompositionalForces(model, forces)
    if ~isempty(forces.W) && numel(forces.W) > 0
        assert(~isempty(model.FacilityModel), ...
            'FacilityModel must be set up before validating driving forces for a compositional problem with wells!');
        T = model.FacilityModel.T;
        p = model.FacilityModel.pressure;
        
        assert(isfield(forces.W, 'components'), ...
            'Wells must have field .components for a compositional model.');
        eos = model.EOSModel;
        if ~isfield(forces.W, 'rhoS')
            act = vertcat(forces.W.status);
            if any(act)
                wellIndices = find(act);
                z = vertcat(forces.W(act).components);
                Z = eos.getMassFraction(z);
                n = size(z, 1);
                [L, x, y, Z_L, Z_V, rhoL, rhoV] = standaloneFlash(repmat(p, n, 1), repmat(T, n, 1), z, eos);

                for i = 1:numel(wellIndices)
                    wNo = wellIndices(i);
                    fprintf('%d of %d\n', i, numel(forces.W));
                    [rho, comp] = getSurfaceParameters(model, forces.W(wNo), rhoL(i), rhoV(i), x(i, :), y(i, :), L(i), Z_L(i), Z_V(i), Z(i, :));
                    forces.W(wNo).compi = comp;
                    forces.W(wNo).rhoS = rho;
                end
            end
        end
    end
end

function [rho, comp] = getSurfaceParameters(model, W, rhoL, rhoV, x, y, L, Z_L, Z_V, Z)
    wat = model.water;
    if false
        % Use flash
        [sL, sV] = eos.computeSaturations(rhoL, rhoV, x, y, L, Z_L, Z_V);
        % compi is a mass-fraction in practice
        L_mass = sL.*rhoL(u)./(sL.*rhoL + sV.*rhoV);
        comp = [L_mass, 1-L_mass];
    else
        % Use the pre-computed definition of light/heavy
        % components to determine "compi"
        isEOS = cellfun(@(x) isa(x, 'EquationOfStateComponent'), model.Components);
        val = cellfun(@(x) x.surfacePhaseMassFractions, model.Components(isEOS), 'UniformOutput', false)';
        val = vertcat(val{:});
        val = val(:, (1+wat):end);
        comp = sum(bsxfun(@times, val, Z'), 1);
    end
    if wat
        assert(~isempty(W.compi), ...
            'W.compi must be present for compositional flow with water phase.');
        sW = W.compi(1);
        comp = [sW, comp.*(1-sW)];
        rho = [model.fluid.rhoWS, rhoL, rhoV];
    else
        rho = [rhoL, rhoV];
    end
end