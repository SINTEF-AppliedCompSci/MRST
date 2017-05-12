function state = addDerivedQuantities(model, state)

    u = model.operators.mech.V_dir; % zeros(G.gridim * G.nodes.num, 1);
    u(~model.operators.mech.isdirdofs) = getProp(model, state, 'xd');
    state = setProp(model, state, 'u', u);
    uu = reshape(u, model.G.griddim, [])';
    state = setProp(model, state, 'uu', uu);
    
    % calculate div, stress, eigen values + + +
    vdiv = model.operators.extra.vdiv * u ./ model.G.cells.volumes;
    state = setProp(model, state, 'vdiv', vdiv);
    
    if(model.G.griddim == 2)
        lin_dim = 3; % dimension of the non trivial linear degree of freedom which
    elseif(model.G.griddim == 3)
        lin_dim = 6;
    else
        error('Wrong dimsension');
    end
    stress = reshape(model.operators.extra.stress * u, lin_dim, [])';
    strain = reshape(model.operators.extra.strain * u, lin_dim, [])';
    
    state = setProp(model, state, 'stress', stress);
    state = setProp(model, state, 'strain', strain);
    
    
end
