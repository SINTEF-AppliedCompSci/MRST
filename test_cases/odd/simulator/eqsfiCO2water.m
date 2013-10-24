function [eqs, info] = eqsfiCO2water(state0, state, dt, Gt, W, s, f, bc, slope, top_temp, temp_grad)

    %% Defining shorthands
    g_cos_t = norm(gravity) * cos(slope.theta); 
    g_sin_t = norm(gravity) * sin(slope.theta);
    H = Gt.cells.H;
    G_cos_t = temp_grad / 1000 * cos(slope.theta);
    G_sin_t = temp_grad / 1000 * sin(slope.theta);
    
    %% current and previous variables
    pI = state.pressure;   pI0 = state0.pressure; % ref. pressure at _interface_
    h  = state.h;          h0  = state0.h;
    
    tI0 = top_temp + h0 * G_cos_t; % initial temperature at _interface_
    
    %% initializing adi-variables
    [pI, h] = initVariablesADI(pI, h);

    tI = top_temp + h * G_cos_t; % current temperature at _interface_ (top_temp a vector here)
    
    %% Setting up densities and correction polynomials

    [Ieta,  INupEta] = etaIntegrals(f.CO2, pI ,  tI, G_cos_t, g_cos_t); 
    Ieta0            = etaIntegrals(f.CO2, pI0, tI0, G_cos_t, g_cos_t);
    
    % densities at interface
    rhoC  = f.CO2.rho(pI, tI);  
    rhoCf = s.faceAvg(rhoC);
    rhoC0 = f.CO2.rho(pI0, tI0);
    rhoW  = f.water.rho;   
    
    %% Establishing caprock and interface gradients
    % (NB: s.grad gives _negative_gradient, so we must use it with an additional '-' sign)
    dpI  = -s.grad(pI);             % grad. of pressure at interface
    dInt = -s.grad(Gt.cells.z + h); % grad. of interface
    
    %% Computing interior fluxes
    v = compute_parall_slopevec(Gt, interiorFaces(Gt), slope.dir);
    iface_modif_term = (g_cos_t .* dInt + g_sin_t .* v);
    CO2_pterm = s.faceAvg(INupEta(-h)) .* (dpI - rhoCf .* iface_modif_term);
    WAT_pterm = dpI - rhoW .* iface_modif_term;
    
    if strcmpi(f.CO2.compressible, 'horizontal')
        CO2_pterm = CO2_pterm - horizontal_correction(g_cos_t, -s.grad(rhoC), s.faceAvg(h));
    end
    
    fluxCO2 = computeFlux(s, f.CO2  , CO2_pterm, rhoC, h,   H);
    fluxWAT = computeFlux(s, f.water, WAT_pterm, rhoW, H-h, H);
    
    %% setting up boundary conditions and well conditions (system right-hand sides)
    
    fix    = find(strcmpi(bc.type, 'flux'));     % ix of flux face indices in bc.face
    pix    = find(strcmpi(bc.type, 'pressure')); % ix of p. face indices in bc.face
    bcells = BFace2Cell(Gt,bc.face);             % cells corresponding to boundary faces

    % establishing right-hand side vectors as ADI variables (this way of
    % doing it might be considered a hack - is there a better way?)
    rhsCO2 = 0 * h; rhsWAT = 0 * h;
    
    % Adding cell in/out fluxes from boundaries with constant-flow conditions 
    % (@@ only supported for water at the moment, so nothing happens with 'rhsCO2')
    rhsWAT = accumarrayADI(bcells(fix), bc.value(fix), rhsWAT); 
    
    % Adding cell in/out fluxes from boundaries with constant-pressure cond.
    bcellsP = bcells(pix);              % boundary cells with pressure-imposed conditions
    bfacesP = bc.face(pix);             % indices of boundary faces
    bpI  = bc.value(pix) + h(bcellsP) * rhoW * g_cos_t;  % boundary press. at interface
    brhoC = f.CO2.rho(bpI, tI(bcellsP)); % density at interface
    bh      = h(bcellsP);               % interface thickness in boundary cells
    bH      = H(bcellsP);               % aquifer thickness at boundary cells
    bdpI = bpI - pI(bcellsP);           % boundary press. grad. at interface, measured
                                        % 'outwards' (NB: factor of 2 omitted, as this
                                        % is already taken into account in
                                        % transmissibilities) 

    % computing parallel gravity components (fixing signs for those faces
    % that have default normals pointing _into_ the cells)
    vI           = compute_parall_slopevec(Gt, bfacesP, slope.dir);
    norm_inv     = find(Gt.faces.neighbors(bfacesP,1) == 0);
    vI(norm_inv) = - vI(norm_inv);

    % compute eta integrals _on_ the pressure boundaries
    [~, bINupEta] = etaIntegrals(f.CO2, bpI, tI(bcellsP), G_cos_t, g_cos_t);
    
    % computing complete expression to go into the divergence term
    % (term in interface gradient skipped here, does it matter?)
    bCO2_pterm = bINupEta(-bh) .* (bdpI - brhoC .* g_sin_t .* vI);
    bWAT_pterm = bdpI - (rhoW * g_sin_t .* vI);
    INupEta_cell_h = INupEta(-h);

    uRhoC    = upstreamVal(rhoC(bcellsP),       brhoC,         bCO2_pterm);
    uINupEta = upstreamVal(INupEta_cell_h(pix), bINupEta(-bh), bCO2_pterm);

    % recomputing bCO2_pterm, ensuring use of upstream values
    bCO2_pterm = uINupEta .* (bdpI - uRhoC .* g_sin_t .* vI);
    
    % Computing boundary fluxes 
    bFluxCO2 = - uRhoC ./ f.CO2.mu   .* (bh ./ bH)    .* s.T_all(bfacesP) .* bCO2_pterm;
    bFluxWAT = -  rhoW ./ f.water.mu .* (bH-bh) ./ bH .* s.T_all(bfacesP) .* bWAT_pterm;
    
    % Adding results to RHS
    rhsCO2 = accumarrayADI(bcellsP, bFluxCO2, rhsCO2);
    rhsWAT = accumarrayADI(bcellsP, bFluxWAT, rhsWAT);

    %% Reporting on CO2 flow across boundary, for logging by higher-level routines

    info.outflowCO2 = valueOf(rhsCO2);       % escaped CO2
    info.outflowWAT = valueOf(rhsWAT);
    info.intPress   = valueOf(pI);           % pressure at interface
    info.topPress   = valueOf(pI + g_cos_t * rhoC .* (-h) .* Ieta(-h)); 
    info.botPress   = valueOf(pI + rhoW .* g_cos_t .* (H-h));
    info.topTemp    = top_temp;
    info.intTemp    = valueOf(tI);
    info.intRho     = valueOf(rhoC);
    info.fluxCO2    = valueOf(fluxCO2);
    info.INupEta    = INupEta(-h);
    if strcmpi(f.CO2.compressible, 'full')
        info.topRho = valueOf(f.CO2.rho(info.topPress, info.topTemp));
    else
        info.topRho = valueOf(rhoC);
    end;

    % Adding fluxes into cells from wells
    [rhsCO2, rhsWAT] = addWellFluxes(rhsCO2, rhsWAT, Gt, W);
    
    %% defining equations

    eqs{1} = (s.pv/dt) .* ((rhoC .* h .* Ieta(-h))  - (rhoC0 .* h0 .* Ieta0(-h0)))  + s.div(fluxCO2) + rhsCO2;    
    eqs{2} = (s.pv/dt) .* rhoW .* (h0 - h) + s.div(fluxWAT) + rhsWAT;

end
%%                   END OF FUNCTION eqsfiCO2water                             

function val = valueOf(var)
    if isa(var, 'ADI')
        val = var.val;
    else
        val = var;
    end
end


%-------------------------------------------------------------------------------
function corr = horizontal_correction(g, grad_rho, h)
%-------------------------------------------------------------------------------
% Compute the correction term needed when considering horizontally-varying density.
    corr = 1/2 * g * grad_rho .* h;
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
    
%-------------------------------------------------------------------------------
function [Cflux, Wflux] = addWellFluxes(Cflux, Wflux, Gt, W)
%-------------------------------------------------------------------------------
    for nw = 1:numel(W)
        w = W(nw);
        assert(strcmpi(w.type, 'rate')); % @@ BHP not yet supported

        % We assume injection is of CO2 (@@injection of water not yet
        % supported)
        % Injection means flow is going _into_ the concerned cells, which is
        % then counted as a negative flux globally, hence the negative sign below.
        Cflux(w.cells) = Cflux(w.cells) - w.val/numel(w.cells);
    end
end
       
%-------------------------------------------------------------------------------
function cells = BFace2Cell(Gt, faces)
%-------------------------------------------------------------------------------
% for each boundary face indexed in the 'faces' vector, identify the unique
% cell that it is a face of.
    cells = sum(Gt.faces.neighbors(faces,:),2); % ix + 0, or 0 + ix
end

%-------------------------------------------------------------------------------
function flux = computeFlux(s, fluid, pterm, rho, h, H)
%-------------------------------------------------------------------------------
    upc = logical(pterm < 0); % identify upstream cells
    % Division by 'H' necessary here, as the transmissibility 's.tI' has been
    % scaled by H 
    flux = -s.faceUpstr(upc, rho .* h ./fluid.mu ./ H) .* s.T .* pterm;
end

%-------------------------------------------------------------------------------
function ifaces = interiorFaces(Gt)
%-------------------------------------------------------------------------------
    % Identify the internal faces (@ cleaner way of doing this?)
    ifaces = find(prod(double(Gt.faces.neighbors), 2) ~= 0);
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

%-------------------------------------------------------------------------------
function ans = upstreamVal(centroid_value, boundary_value, grad)
%-------------------------------------------------------------------------------    

    ans = centroid_value;

    % If gradient is positive, pressure is increasing towards boundary, which
    % means that the flow is directed into the cell.  In that case, the
    % upstream value is the value on the boundary.  Otherwise, the opposite
    % is the case.
    pos_grad = grad > 0;    
    ans(pos_grad) = boundary_value(pos_grad);
end

% %-------------------------------------------------------------------------------
% function res = compute_pseudo_g(g, rho, dtop)
% %-------------------------------------------------------------------------------
 
%     res =  g * dtop .* rho;
% end

