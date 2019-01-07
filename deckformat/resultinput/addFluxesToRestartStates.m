function states = addFluxesToRestartStates(states, G, T, fluid, deck)
gravity reset on
% get model with empty rock for rel-perm computations
model = selectModelFromDeck(G, [], fluid, deck);
model.fluid = fluid;
% setup only the operators we need
model.operators = setupOperatorsTPFA(G, [], 'trans', sparse(numel(T),1), 'porv', sparse(G.cells.num,1));
gdz = model.getGravityGradient();

if ~iscell(states)
    states = {states};
end
for k =1:numel(states)
    % relperm
    if model.water && model.oil  && model.gas
        [kr{1}, kr{2}, kr{3}] = model.evaluateRelPerm({states{k}.s(:,1), states{k}.s(:,2), states{k}.s(:,3)});
    else
        [kr{1}, kr{2}] = model.evaluateRelPerm({states{k}.s(:,1), states{k}.s(:,2)});
    end
    % fluxes
    v = cell(1,3);
    phNum = 0;
    if model.water
        phNum = phNum +1;
        v{phNum} = getFluxAndPropsWater_BO(model, states{k}.pressure, states{k}.s(:,phNum), kr{phNum}, T, gdz);
    end
    if model.oil
        phNum = phNum +1;
        if model.disgas
            v{phNum} = getFluxAndPropsOil_BO(model, states{k}.pressure, states{k}.s(:,phNum), kr{phNum}, T, gdz, ...
                                         states{k}.rs,  states{k}.s(:,phNum+1)>0);
        else
            v{phNum} = getFluxAndPropsOil_BO(model, states{k}.pressure, states{k}.s(:,phNum), kr{phNum}, T, gdz);
        end
    end
    if model.gas
        phNum = phNum +1;
        if model.vapoil
            v{phNum} = getFluxAndPropsGas_BO(model, states{k}.pressure, states{k}.s(:,phNum), kr{phNum}, T, gdz, ...
                                         states{k}.rv,  states{k}.s(:,phNum-1)>0);
        else
            v{phNum} = getFluxAndPropsGas_BO(model, states{k}.pressure, states{k}.s(:,phNum), kr{phNum}, T, gdz);
        end
    end
    states{k} = model.storeFluxes(states{k}, v{1}, v{2}, v{3});
end
end