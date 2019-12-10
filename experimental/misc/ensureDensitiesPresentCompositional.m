function state = ensureDensitiesPresentCompositional(model, state)
model.FacilityModel = FacilityModel(model);
state = model.computeFlash(state, inf);
state.wellSol = [];
model.extraStateOutput = true;
[~, state] = model.getEquations(state, state, 1, struct('W', [], 'src', [], 'bc', []));

end