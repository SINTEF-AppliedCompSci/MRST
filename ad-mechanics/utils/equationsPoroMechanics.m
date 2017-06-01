function [eqs, names, types, state] = equationsPoroMechanics(x, model, fluidp)
    
    G = model.G;
    s = model.operators.mech;
    alpha = model.rock.alpha;

    eqs{1} = s.A * x - s.gradP * (alpha .* fluidp) - s.rhs;

    fac =  1 / (1e6 * mean(G.cells.volumes));
    eqs{1} = eqs{1} * fac;
    names = {'disp'};
    types = {'disp_dofs'};

end