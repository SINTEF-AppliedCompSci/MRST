function state = assignDofFromState(disc, state)

    nDof = disc.basis.nDof;
    G    = disc.G;
    sdof = zeros(G.cells.num*nDof, size(state.s,2));
    sdof(1:nDof:G.cells.num*nDof,:) = state.s;
    
    state.sdof = sdof;
    
end