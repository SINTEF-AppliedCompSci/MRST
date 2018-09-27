function [flag_v, flag_G, upCells_v, upCells_G, s_v, s_G] = getSaturationUpwindDG(disc, faces, x, T, vT, g, mob, sdof, state)
    % Explicit calculation of upstream cells (Bernier & Jaffre)
    % for each quadrature point x on each face in faces.
    %
    % PARAMETERS:
    %   faces - Faces we want upstream cells for
    %   x     - Upstram cells are computed for face(ix) at x(ix,:)
    %   sdof  - Saturation degrees of freedom
    %   state - For dofPos
    %   T, vT, g, mob - Inputs to getSaturationUpwind
    %
    % RETURNS:
    %   flag_v, flag_G       - Upstream flags for viscous and gravity
    %   upCells_v, upCells_G - Upstream cells for viscous and gravity
    %   s_v, sG              - Corresponding saturations*
    %
    %   Returned from function since since they must anyway be
    %   calculated later on

    G = disc.G;

    % Mapping from all faces to internal connections
    all2int = zeros(G.faces.num,1);
    all2int(disc.internalConn) = 1:nnz(disc.internalConn);
    ix      = all2int(faces);

    % Extract cells, trnasmissibilities velocities, total fluxed
    % and gravity terms
    cL = disc.N(ix,1);
    cR = disc.N(ix,2);
    T  = T(ix);
    vT = vT(faces);
    g  = cellfun(@(g) g(faces), g, 'unif', false);

    % Transform to referece coordinates
    [xL, ~, ~] = disc.transformCoords(x, cL);
    [xR, ~, ~] = disc.transformCoords(x, cR);

    % Evaluate saturations and mobilities on each side of each face
    % at each quadrature point
    s = disc.evaluateSaturation([xL; xR], [cL; cR], sdof, state);
    sL = s(1:numel(cL)); sR = s(numel(cL)+1:end);
    mob{1} = mob{1}([sL; sR], [cL; cR]);
    mob{2} = mob{2}(1-[sL; sR], [cL; cR]);

    % Make fake faceUpstr function
    N = [1:numel(ix); numel(ix)+1:2*numel(ix)]';
    upw = @(flag, x)faceUpstr(flag, x, N, [size(N,1), max(max(N))]);

    % Use standard MRST function to compute upstraeam flags
    [flag_v, flag_G] = getSaturationUpwind(disc.upwindType, state, g, vT, T, mob, upw);

    % For each phase, assign upstram cell and corresponding
    % saturation for each quadrature point of each face
    nPh = numel(g);
    [upCells_v, upCells_G] = deal(repmat({cR}, 1, nPh));
    [s_v, s_G] = deal(cell(1, nPh));
    [s_v{:},s_G{:}] = deal(sR);
    for phNo = 1:nPh
        % Viscous upstream cell
        upCells_v{phNo}(flag_v(:,phNo)) = cL(flag_v(:,phNo));
        s_v{phNo}(flag_v(:,phNo)) = sL(flag_v(:,phNo));
        % Gravity upstram cell
        upCells_G{phNo}(flag_G(:,phNo)) = cL(flag_G(:,phNo));
        s_G{phNo}(flag_G(:,phNo)) = sL(flag_G(:,phNo));
    end

end