function assembly = assembleBiot(G, props, drivingforces, eta, tbls, mappings, varargin)
    
    mechprops = props.mechprops;
    K = props.fluidprops;
    
    % Assemble mechanic problem
    mechforces = drivingforces.mechanics;
    loadstruct = mechforces.loadstruct;
    mechassembly = assembleMPSA(G, mechprop, loadstruct, eta, tbls, mappings)
    
    % Assemble fluid problem
    fluidforces = drivingforces.fluid;
    bcstruct = fluid.bcstruct;
    src      = fluid.src;
    fluidassembly = assembleMPFA(G, K, bcstruct, src, eta, tbls, mappings, varargin)

    % Compute coupling terms
    
    % Recover matrices from mechanic assembly
    mechmat = mechassembly.matrices;
    invA11 = mechmat.invA11;
    A12 = mechmat.A12;
    A21 = mechmat.A21;
    A15 = -mechmat.D;
    A51 = -A15';
    
    % Recover matrices from mechanic assembly
    fluidmat = fluidassembly.matrices;
    invA33 = fluidmat.invA11;
    A34 = fluidmat.A12;
    A43 = fluidmat.A43;
    A36 = -fluidmat.D;
    A63 = -A36';
    
    % Assemble the coupling terms (A14 and A41)
    
    % 1. row : Momentum equation
    A21invA11 = A21*invA11;
    B11 = -A21invA11*A12 +  A22;
    B12 = -A21invA11*A14;
    B13 = -A21invA11*A15;

    % 2. row : Fluid mass conservation
    A41invA11 = A41*invA11;
    A43invA33 = A43*invA33;
    B21 = -A41invA11*A12;
    B22 = -A41invA11*A14 - A43invA33*A34 + A44;
    B23 = -A41invA11*A15;
    B24 = -A43invA33*A36;

    % 3. row : Mechanic BC
    A51invA11 = A51*invA11;
    B31 = -A51invA11*A12;
    B32 = -A51invA11*A14;
    B33 = -A51invA11*A15;

    % 4. row : Fluid BC
    A63invA33 = A63*invA33;
    B42 = -A63invA33*A34;
    B42 = -A63invA33*A36;


end

