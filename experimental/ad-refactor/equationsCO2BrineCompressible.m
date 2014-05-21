function [problem, state] = equationsCO2BrineCompressible(state0, state, dt, ...
                                                          Gt, drivingForces, s, ...
                                                          co2fluid, brinefluid, ...
                                                          T_ref, T_ref_depth, T_grad, ...
                                                          slope, slopedir, ...
                                                          varargin)
    opt = struct('reverseMode' , false, 'resOnly', false, 'iteration', -1);
    opt = merge_options(opt, varargin{:});
    
    %% Ensure that the state's wellSol struct is present if and only if there
    %% are wells present
    if isempty(drivingForces.Wells) && isfield(state, 'wellSol')
        % no wells present - remove wellSol structure
        state.wellSol.qGs = [];
        state.wellSol.bhp = [];
    elseif ~isempty(drivingForces.Wells) && ~isfield(state, 'wellSol')
         % wells present and no wellSol structure initialized.  Do it.
         state.wellSol = struct('qGs', drivingForces.Wells.val, ...
                                'bhp', state.pressure(vertcat(drivingForces.Wells.wc)));
    end
        
    %% Precomputed values useful for later
    g_cos_t = norm(gravity) * cos(slope);
    g_sin_t = norm(gravity) * sin(slope);
    H = Gt.cells.H;
    
    %% Extract current and prevous variables, initialize as ADI if required
    [pI,  h,  q,  bhp] = extract_vars(state, ~opt.resOnly && ~opt.reverseMode);
    [pI0, h0, ~,  ~  ] = extract_vars(state0,~opt.resOnly &&  opt.reverseMode);
    
    %% Computing values at co2-brine interface
    % temperature
    tI  = imposed_iface_temp(Gt, slope, slopedir, h,  T_ref, T_ref_depth, T_grad);
    tI0 = imposed_iface_temp(Gt, slope, slopedir, h0, T_ref, T_ref_depth, T_grad);
    % density and viscosity at interface conditions
    rhoC  = co2fluid.rho(pI,  tI)  ;   rhoB  = brinefluid.rho(pI,  tI);
    rhoC0 = co2fluid.rho(pI0, tI0) ;   rhoB0 = brinefluid.rho(pI0, tI0);
    muC   = co2fluid.mu(pI, tI)    ;   muB   = brinefluid.mu(pI, tI);
    
    %% Setting up density-corrective integrals
    [IetaCO2, INupEtaCO2] = co2fluid.h_integrals(pI, tI);
    [Ieta0CO2, ~]         = co2fluid.h_integrals(pI0, tI0);
    [IetaBri, INupEtaBri] = brinefluid.h_integrals(pI, tI);
    [Ieta0Bri, ~]         = brinefluid.h_integrals(pI0, tI0);
    
    %% Computing interior fluxes
    
    % NB: s.grad gives _negative_gradient, hence the additional '-' sign) 
    dpI  = -s.grad(pI);              % gradient of interface pressure
    dInt = -s.grad(Gt.cells.z + h);  % gradient of interface position
    
    % combined gravity/interface slope term
    v = compute_parall_slopevec(Gt, interiorFaces(Gt), slopedir);
    mod_term = (g_cos_t .* dInt + g_sin_t .* v);
    
    intFluxCO2 = computeFlux(s.faceAvg, s.faceUpstr, INupEtaCO2, s.T, -1, ...
                             dpI, mod_term, -s.grad(rhoC), rhoC, muC, h,   H, g_cos_t); 
    intFluxBri = computeFlux(s.faceAvg, s.faceUpstr, INupEtaBri, s.T,  1, ...
                             dpI, mod_term, -s.grad(rhoB), rhoB, muB, H-h, H, g_cos_t);
    
    %% Computing boundary fluxes, if any
    [bc_cell, bc_fC, bc_fB] = ...
        BCFluxes(Gt, s, drivingForces.bc, pI, h, slope, slopedir, ...
                 rhoC, muC, INupEtaCO2, rhoB, muB, INupEtaBri);
    
    %% Handling wells (setting up equations, rates)
    [eqs(3:4), wc, qw, types(3:4), names(3:4)] = ...
        setupWellEquations(drivingForces.Wells, pI, q, bhp, rhoC, muC);
    
    %% Setting up equation systems
    rhoCIeta = rhoC; rhoCIeta0 = rhoC0; rhoBIeta = rhoB; rhoBIeta0 = rhoB0;
    if ~isempty(IetaCO2) 
        rhoCIeta  = rhoCIeta  .* IetaCO2(-h);  
        rhoCIeta0 = rhoCIeta0 .* Ieta0CO2(-h0);  
    end
    if ~isempty(IetaBri)
        rhoBIeta  = rhoBIeta  .* IetaBri(H-h);
        rhoBIeta0 = rhoBIeta0 .* Ieta0Bri(H-h);
    end
    
    eqs{1} = contEquation(s, dt, h,   h0,   rhoCIeta, rhoCIeta0, intFluxCO2, bc_cell, bc_fC, wc, qw);
    eqs{2} = contEquation(s, dt, H-h, H-h0, rhoBIeta, rhoBIeta0, intFluxBri, bc_cell, bc_fB, [], []);
    
    types(1:2) = {'cell', 'cell'};
    names(1:2) = {'CO2', 'brine'};

    % Rebalancing the conservation equations
    for i = [1:2]   eqs{i} = eqs{i} * dt/year;   end %#ok

    primary = {'pressure' ,'height', 'q', 'bhp'};
    problem = linearProblem(eqs, types, names, primary, state);
    
    %% Adding extra information to state object (not used in equations, but
    %% useful for later analysis)
    state.extras = addedInfoForAnalysis(tI, rhoC, intFluxCO2, intFluxBri);

    % if ~isempty(IetaCO2)    krull1 = IetaCO2(-h);    else krull1 = ones(100,1); end;
    % if ~isempty(INupEtaCO2) krull2 = INupEtaCO2(-h); else krull2 = ones(100,1); end;
    %state.extras.h50 = [double(state.h(50)), double(krull1(50)), double(krull2(50))]; @@ REMOVE
end

% ----------------------------------------------------------------------------
function extras = addedInfoForAnalysis(tempIface, rhoI, fluxCO2, fluxBrine)
% ----------------------------------------------------------------------------
    
    extras.tI = double(tempIface);
    extras.rhoI = rhoI;
    extras.fluxCO2 = double(fluxCO2);
    extras.fluxBrine = double(fluxBrine);
end
    

% ----------------------------------------------------------------------------
function [eqs, wcells, wrates, tp, nm] = setupWellEquations(W, p, q, bhp, rhoC, muC)
% ----------------------------------------------------------------------------
    if isempty(W)  
        eqs{1} = 0*bhp; 
        eqs{2} = 0*bhp; 
        wcells = [];
        wrates = [];
        tp = {'none', 'none'};
        nm = {'empty', 'empty'};
        return;
    end
    
    % at the moment, only injectors of CO2 are allowed
    assert(all(cellfun(@(x)strcmp('rate', x), {W.type}')));
    
    % Use getWellStuffOG to get the (limited) information we need
    [Tw, ~, Rw, wcells, perf2well] = getWellStuffOG(W);
    
    wrates = rhoC(wcells) .* (-Tw) ./ muC(wcells) .* (p(wcells) - bhp(perf2well));
    eqs{1} = -Rw' * wrates + q;  % @@ Why has the other code a +zeroW here?
    eqs{2} = handleBC(W, bhp, [], [], q);
    tp = {'perf', 'well'};
    nm = {'gasWells', 'closureWells'};
end

% ----------------------------------------------------------------------------
function eq = contEquation(s, dt, h, h0, rho, rho0, intflux, bcells, bflux, wcells, wrate)
% Function to construct the continuity equation.
% ----------------------------------------------------------------------------
    eq = (s.pv/dt) .* ((rho .*h) - (rho0 .* h0)) + s.div(intflux);
    
    eq(wcells) = eq(wcells) - wrate; % correct for wells 

    % @@ NB: The code line below is not correct for cells with more than one
    % boundary face with pressure condition.  'accumarray' should really be
    % used, but the ADI framework currently does not support this function,
    % and a manual loop ('accumarrayADI') is way to slow (practically doubles
    % the whole runtime of the simulation), so we neglect the error for now.
    eq(bcells) = eq(bcells) + bflux; % correct for boundary conditions
end


%-------------------------------------------------------------------------------
function [cells, fluxC, fluxB] = BCFluxes(Gt, s, bc, p, h, slope, slopedir, ...
                                          rhoC, muC, INEfunC, rhoB, muB, INEfunB)
%-------------------------------------------------------------------------------
    if isempty(bc);
        % All boundary conditions are no-flow.  Nothing to do here.
        cells = []; fluxC = []; fluxB = []; return;
    end
    assert(all(strcmp(bc.type, 'pressure'))); % only support pressure type for now

    
    Tbc   = s.T_all(bc.face);
    cells = sum(Gt.faces.neighbors(bc.face, :), 2);

    %assert(numel(unique(cells)) == numel(cells)); % only support one bc per 

    % prepare boundary-cell-specific values
    bdp = bc.value - p(cells);
    bh  = h(cells)     ;  bH  = Gt.cells.H(cells);
    brC = rhoC(cells)  ;  brB = rhoB(cells);
    bmC = muC(cells)   ;  bmB = muB(cells);

    % Computing the slope across boundary faces, and correcting for faces
    % whose face normals don't point 'out'. 
    bv         = compute_parall_slopevec(Gt, bc.face, slopedir);
    inv_ix     = find(Gt.faces.neighbors(bc.face,1) == 0);
    bv(inv_ix) = -bv(inv_ix);

    % Computing the modification term (ignoring caprock shape at boundary,
    % since we cannot really know that...)
    mterm = norm(gravity) * sin(slope) * bv;

    % Manually setting up the relevant discretization operators for the boundary
    avg   = @(x) x;     % x will already refer to boundary faces here, so nothing to do.
    ustr  = @(i, x) x; % ditto
    bINEfunC = make_restrict_fun(INEfunC, cells, rhoC); % @@ Suboptimal solution.  Candidate
    bINEfunB = make_restrict_fun(INEfunB, cells, rhoB); %    for optimization if necessary.
    gct   = norm(gravity) * cos(slope);    
    fluxC = computeFlux(avg, ustr, bINEfunC, Tbc, -1, bdp, mterm, 0, brC, bmC, bh, bH, gct);
    fluxB = computeFlux(avg, ustr, bINEfunB, Tbc,  1, bdp, mterm, 0, brB, bmB, bH-bh, bH, gct);

end

% ----------------------------------------------------------------------------
function res = make_restrict_fun(fun, ix, model)
% @@ A really hacky function.  It would be better to change the calling
% code so that the below isn't necessary.
% ----------------------------------------------------------------------------
    if isempty(fun)
        res = [];
        return;
    end
    
    function v = f(x)
        tmp = model*0; % We just need an empty (possible ADI) variable of
                       % the size expected by the function 'fun'
        tmp(ix) = x;   % We only need the evaluations of 'fun' for the        
        v = fun(tmp);  % indices 'ix'.  Give these indices a non-dummy value, 
        v = v(ix);     % compute the function, and extract only the resulting 
    end                % components.                                          
        
    res = @f;
end


%-------------------------------------------------------------------------------
function flux = computeFlux(faceAvgFun, UpstrFun, INupEtaFun, ...
                            trans, dir, dpI, modif_term, drho, rho, mu, h, H, gct)
%-------------------------------------------------------------------------------  
    
    pterm = dpI - faceAvgFun(rho) .* modif_term;
    
    if ~isempty(INupEtaFun)
        % Apply full compressibility correction in depth
        pterm = faceAvgFun(INupEtaFun(dir*h)) .* pterm;
    else
        % constant density in depth - apply semi-compressible correction factor instead
        pterm = pterm + horizontal_correction(gct, drho, faceAvgFun(dir*h));
    end
    
    upc  = logical(pterm < 0);

    % Division by 'H' necessary here, as the transmissibility 's.tI' has been
    % scaled by H 
    flux = -UpstrFun(upc, rho .* h ./ mu ./ H) .* trans .* pterm;
    
end

%-------------------------------------------------------------------------------
function corr = horizontal_correction(g, grad_rho, h)
%-------------------------------------------------------------------------------
% Compute the correction term needed when considering horizontally-varying density.
    corr = 1/2 * g * grad_rho .* h;
end

%-------------------------------------------------------------------------------
function slope2D = compute_parall_slopevec(Gt, face_ix, slopedir)
%-------------------------------------------------------------------------------
    
    % 'slope2D' shall represent the component of the true 'down' vector that
    % is parallel to the slope.  The argument 'slopedir' is however pointing
    % in the opposite direction (the direction of gravity-driven flow), so we
    % must use a negative sign to invert it.
    v2D_unit = - slopedir ./ norm(slopedir);
        
    %   To find the component that is perpendicular to
    %   a given face, we must take the scalar product with the normal to 
    % face.  However, the normals stored in Gt are not unitary, but scaled
    % according to the areas of the respective faces.  We thus need to divide
    % by face areas in order to get the unitary normals.  
    %   On the other hand, in order to remain consistent with the way spatial
    % gradients are treated in in the equation above, we should multiply by
    % the distance between the centroids of the two cells (or boundaries) involved.
    %   Thus, we here compute a set of internal normals that are divided by
    % face area and multiplied by distance between cells (or boundary faces)
    neigh_cells = Gt.faces.neighbors(face_ix,:);

    lcenter = zeros(numel(face_ix), 2);
    rcenter = zeros(numel(face_ix), 2);
    
    b_at_left  = (neigh_cells(:,1) == 0);
    b_at_right = (neigh_cells(:,2) == 0);

    lcenter( b_at_left,:)  = Gt.faces.centroids(face_ix(b_at_left),          :);
    rcenter( b_at_right,:) = Gt.faces.centroids(face_ix(b_at_right),         :);
    lcenter(~b_at_left,:)  = Gt.cells.centroids(neigh_cells(~b_at_left,  1), :);
    rcenter(~b_at_right,:) = Gt.cells.centroids(neigh_cells(~b_at_right, 2), :);
    
    neigh_dists = row_norms(lcenter - rcenter);
    
    normals = Gt.faces.normals(face_ix,:);
    scale = neigh_dists ./ row_norms(normals);
    normals = diag(scale) * normals;
    
    slope2D = normals * v2D_unit';
    
end

% ----------------------------------------------------------------------------
function tI = imposed_iface_temp(Gt, slope, slopedir, h, T_ref, T_ref_depth, T_grad)
% ----------------------------------------------------------------------------
    depth = computeRealDepth(Gt, slope, slopedir, h);
    tI    = T_ref + (depth - T_ref_depth) * T_grad/1000;
end

% ----------------------------------------------------------------------------
function [p, h, q, bhp] = extract_vars(st, initADI)
% ----------------------------------------------------------------------------
    
    [p, h] = deal(st.pressure, st.h);
    if ~isempty(st.wellSol)
        [q, bhp] = deal(vertcat(st.wellSol.qGs), ...
                        vertcat(st.wellSol.bhp));
    else
        q = 0; bhp = 0;
    end
    

    if initADI
        [p, h, q, bhp] = initVariablesADI(p, h, q, bhp);
    end
end

%-------------------------------------------------------------------------------
function ifaces = interiorFaces(Gt)
%-------------------------------------------------------------------------------
    % Identify the internal faces (@ cleaner way of doing this?)
    ifaces = find(prod(double(Gt.faces.neighbors), 2) ~= 0);
end


%-------------------------------------------------------------------------------
function target = accumarrayADI(subs, val, target)
%-------------------------------------------------------------------------------
% @@ This can be considered a 'slow hack', to get around the problem that
% accumarray does not work with ADI variables.  Is there a faster way of
% doing this, rather than using an explicit loop?
    for i = 1:numel(subs)
        ix  = subs(i);
        inc =  val(i); 
        target(ix) = target(ix) + inc;
    end
end    
