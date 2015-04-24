function [problem, state] = equationsCO2BrineCompressible(state0, state, dt, ...
                                                          Gt, drivingForces, s, ...
                                                          fluid, ...
                                                          T_ref, T_ref_depth, T_grad, ...
                                                          slope, slopedir, ...
                                                          varargin)
   % NB, 'pressure_lock' options only intended for use on systems that would
   % otherwise have indeterminate pressure, i.e. incompressible systems with
   % no pressure boundary conditions.
    opt = struct('reverseMode' , false, 'resOnly', false, 'iteration', -1, ...
                 'pressure_lock', []);
    opt = merge_options(opt, varargin{:});
    
    
    %% Ensure that the state's wellSol struct is present if and only if there
    %% are wells present
    if isempty(drivingForces.Wells) && isfield(state, 'wellSol') && ~isempty(state.wellSol)
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
    rhoC  = fluid.gas.rho(pI,  tI)  ;   rhoB  = fluid.wat.rho(pI,  tI);
    rhoC0 = fluid.gas.rho(pI0, tI0) ;   rhoB0 = fluid.wat.rho(pI0, tI0);
    muC   = fluid.gas.mu(pI, tI)    ;   muB   = fluid.wat.mu(pI, tI);
    
    %% Setting up density-corrective integrals
    [IetaCO2, INupEtaCO2, ~] = fluid.gas.h_integrals(pI, tI);
    [Ieta0CO2, ~, ~]         = fluid.gas.h_integrals(pI0, tI0);
    [IetaBri, INupEtaBri, ~] = fluid.wat.h_integrals(pI, tI);
    [Ieta0Bri, ~, ~]         = fluid.wat.h_integrals(pI0, tI0);
    
    %% Computing interior fluxes
    dpI  = s.Grad(pI);              % gradient of interface pressure
    dInt = s.Grad(Gt.cells.z + h);  % gradient of interface position

    % Override periodic boundary condition faces, if any
    if isfield(drivingForces, 'bcp') && ~isempty(drivingForces.bcp)
       iface_ixs = interiorFaces(Gt);
       reindex(iface_ixs) = 1:numel(iface_ixs);
       bcp = drivingForces.bcp; % for convenience
       dpI(reindex(bcp.face)) = dpI(reindex(bcp.face)) + bcp.value .* bcp.sign;
    end

    % combined gravity/interface slope term
    mod_term = g_cos_t .* dInt;
    if norm(g_sin_t) > 0
       v        = compute_parall_slopevec(Gt, interiorFaces(Gt), slopedir);
       mod_term = mod_term + g_sin_t .* v;
    end
    
    intFluxCO2 = computeFlux(s.faceAvg, s.faceUpstr, INupEtaCO2, fluid.gas.kr, s.T, -1, ...
                             dpI, mod_term, s.Grad(rhoC), rhoC, muC, h, H, g_cos_t); 
    intFluxBri = computeFlux(s.faceAvg, s.faceUpstr, INupEtaBri, fluid.wat.kr, s.T,  1, ...
                             dpI, mod_term, s.Grad(rhoB), rhoB, muB, H-h, H, g_cos_t);
    
    %% Computing boundary fluxes, if any
    [bc_cell, bc_fC, bc_fB] = ...
        BCFluxes(Gt, fluid, s, drivingForces.bc, pI, h, slope, slopedir, ...
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

    if ~isempty(opt.pressure_lock)
       % Replacing one of the equations of the presumably degenerate system
       % with an equation linking the pressure of a cell to the requested pressure
       eqs{1}(1) = pI(1) - opt.pressure_lock;
    end
    
    % Rebalancing the conservation equations
    for i = [1:2]   eqs{i} = eqs{i} * dt/year;   end %#ok

    primary = {'pressure' ,'height', 'qGs', 'bhp'};
    problem = LinearizedProblem(eqs, types, names, primary, state);
    problem.iterationNo = opt.iteration; % @ should ideally be unnecessary
    
    %% Adding extra information to state object (not used in equations, but
    %% useful for later analysis)
    state.extras = addedInfoForAnalysis(tI, rhoC, intFluxCO2, intFluxBri);
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
    
    % Use getWellStuffWG to get the (limited) information we need
    [Tw, ~, Rw, wcells, perf2well] = getWellStuffWG(W);
    
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
    eq = (s.pv/dt) .* ((rho .*h) - (rho0 .* h0)) + (s.Div(intflux));
    
    eq(wcells) = eq(wcells) - wrate; % correct for wells 

    % @@ NB: The code line below is not correct for cells with more than one
    % boundary face with pressure condition.  'accumarray' should really be
    % used, but the ADI framework currently does not support this function,
    % and a manual loop ('accumarrayADI') is way to slow (practically doubles
    % the whole runtime of the simulation), so we neglect the error for now.
    eq(bcells) = eq(bcells) + bflux; % correct for boundary conditions
end


%-------------------------------------------------------------------------------
function [cells, fluxC, fluxB] = BCFluxes(Gt, fluid, s, bc, p, h, slope, slopedir, ...
                                          rhoC, muC, INEfunC, rhoB, muB, INEfunB)
%-------------------------------------------------------------------------------
    if isempty(bc);
        % All boundary conditions are no-flow.  Nothing to do here.
        cells = []; fluxC = []; fluxB = []; return;
    end
    assert(all(strcmp(bc.type, 'pressure'))); % only support pressure type for now

    
    Tbc   = s.T_all(bc.face);
    cells = sum(Gt.faces.neighbors(bc.face, :), 2);
    nc    = numel(cells);

    %assert(numel(unique(cells)) == numel(cells)); % only support one bc per 
    % prepare boundary-cell-specific values (first: face values, then: cell values)
    bh  = repmat(h(cells),    2, 1)  ;  bH  = repmat(Gt.cells.H(cells), 2, 1);
    brC = repmat(rhoC(cells), 2, 1)  ;  brB = repmat(rhoB(cells),       2, 1);
    bmC = repmat(muC(cells),  2, 1)  ;  bmB = repmat(muB(cells),        2, 1);
    if ~isempty(bc.sat)
       % Adjusting h-value at boundary according to prescribed saturation there.
       bh(1:nc) = bc.sat .* bH(1:nc);
    end
    
    bp  = bc.value + brB(1:nc) * norm(gravity) .* bh(1:nc) * cos(slope);
    bdp = bp - p(cells);

    % Computing the slope across boundary faces, and correcting for faces
    % whose face normals don't point 'out'. 
    if norm(slope) > 0
       bv = compute_parall_slopevec(Gt, bc.face, slopedir);
    else
       bv = zeros(numel(bc.face), 1);
    end
    
    inv_ix     = find(Gt.faces.neighbors(bc.face,1) == 0);
    bv(inv_ix) = -bv(inv_ix);

    % Computing the modification term (ignoring caprock shape at boundary,
    % since we cannot really know that...)
    mterm = (bh(1:nc) - bh(nc+1:2*nc)) + norm(gravity) * sin(slope) * bv;

    % Manually setting up the relevant discretization operators for the boundary
    avg      = @(x)    bcond_average(x, nc);
    ustr     = @(i, x) bcond_ustr(i, x, nc);
    bINEfunC = make_restrict_fun(INEfunC, cells, rhoC); % @@ Suboptimal solution.  Candidate
    bINEfunB = make_restrict_fun(INEfunB, cells, rhoB); %    for optimization if necessary.
    bkrG     = make_restrict_fun(fluid.gas.kr, cells, rhoC);
    bkrW     = make_restrict_fun(fluid.wat.kr, cells, rhoB);
    gct      = norm(gravity) * cos(slope);
    
    % Computing fluxes
    fluxC = computeFlux(avg, ustr, bINEfunC, bkrG, Tbc, -1, bdp, mterm, 0, brC, bmC, bh, bH, gct);
    fluxB = computeFlux(avg, ustr, bINEfunB, bkrW, Tbc,  1, bdp, mterm, 0, brB, bmB, bH-bh, bH, gct);
end

% ----------------------------------------------------------------------------
function res = bcond_ustr(i, x, numel)
% ----------------------------------------------------------------------------
%   res = x(numel+1:2*numel);
   res     = x(1:numel);         % face values
   cvals   = x(numel+1:2*numel); % make separate vector of cell values
   res(i) = cvals(i);            % insert cell values into result vector where required
end

% ----------------------------------------------------------------------------
function res = bcond_average(x, numel)
% ----------------------------------------------------------------------------
%res = x(numel+1: 2*numel);
res = 0.5 * (x(1:numel) + x(numel+1:2*numel));
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
    
    % function v = f(x)
    %     tmp = model*0; % We just need an empty (possible ADI) variable of
    %                    % the size expected by the function 'fun'
    %     tmp(ix) = x;   % We only need the evaluations of 'fun' for the        
    %     v = fun(tmp);  % indices 'ix'.  Give these indices a non-dummy value, 
    %     v = v(ix);     % compute the function, and extract only the resulting 
    % end                % components.                                          
    
    N = numel(ix);    
    function v = f(x)
       tmp_faces = model * 0;
       tmp_cells = model * 0;
       tmp_faces(ix) = x(1:N);
       tmp_cells(ix) = x(N+1:2*N);
       v_faces = fun(tmp_faces);
       v_cells = fun(tmp_cells);
       v = [v_faces(ix); v_cells(ix)];
    end
    
    
    res = @f;
end


%-------------------------------------------------------------------------------
function flux = computeFlux(faceAvgFun, UpstrFun, INupEtaFun, relPermFun,...
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
    flux = -UpstrFun(upc, rho .* relPermFun(h./H) ./ mu) .* trans .* pterm;
    
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
    if isfield(st, 'wellSol') && ~isempty(st.wellSol)
        [q, bhp] = deal(vertcat(st.wellSol.qGs), ...
                        vertcat(st.wellSol.bhp));
    else
        q = []; bhp = [];
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

% ----------------------------------------------------------------------------
function res = ifelse(cond, yes, no)
% ----------------------------------------------------------------------------
   if (cond) res = yes; else yes = no; end
end
