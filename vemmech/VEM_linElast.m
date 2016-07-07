function [uu, extra] = VEM_linElast(G, C, el_bc, load, varargin)
%
%
% SYNOPSIS:
%   function [uu, extra] = VEM_linElast(G, C, el_bc, load, varargin)
%
% DESCRIPTION: solve linear elastisity problem using the Virtual Element method. Main
% mechanics solver. Without face degree of freedom. Meant to be run for any
% dimension.
%
% PARAMETERS:
%   G        - Grid structure as described by grid_structure, which also
%              has the mappings from G = createAugmentedGrid(G).
%              Some options need G = computeGeometryCalc(G).
%   C        - Elasticity tensor
%   el_bc    - boundary condition of type struct('disp_bc', [], 'force_bc', [])
%   load     - loading term
%   varargin -
%
% RETURNS:
%   uu    - Displacement field
%   extra - Extra outputs
%
% EXAMPLE:
%
% COMMENTS: The assembly code is for general elastisty tensor in voit's notation:
% http://en.wikipedia.org/wiki/Voigt_notation
%
% SEE ALSO:


    opt = struct('type'         , 'mrst'                          , ...
                 'linsolve'     , @mldivide                       , ...
                 'blocksize'    , 30000 - (G.griddim - 2) * 20000 , ...
                 'add_operators', true                            , ...
                 'force_method' , 'dual_grad_type'                , ...
                 'alpha_scaling', 1                               , ...
                 'S'            , []                              , ...
                 'experimental_scaling', false                    , ...
                 'pressure'     , []);
    opt = merge_options(opt, varargin{:});
    opt.add_operators = opt.add_operators && nargout>1;
    if(opt.add_operators)
        [S, extra] = VEM_mrst_vec(G, C, ...
                                  'blocksize'            , opt.blocksize, ...
                                  'alpha_scaling'        , opt.alpha_scaling, ...
                                  'S'                    , opt.S, ...
                                  'experimental_scaling' , opt.experimental_scaling);
    else
        S = VEM_mrst_vec(G, C, ...
                         'blocksize'           , opt.blocksize, ...
                         'alpha_scaling'       , opt.alpha_scaling, ...
                         'S'                   , opt.S, ...
                         'experimental_scaling', opt.experimental_scaling);
    end

    % Recalculate "weights" (all are calculated in the assembly, they could in fact only be calculated
    % once).
    if (G.griddim == 2)
        qf_all = G.faces.areas / 2;
        [qc_all, qcvol] = calculateQF(G);
    else
        [qc_all, qf_all, qcvol] = calculateQC(G);
    end

    % Apply Diriclet boundary conditions
    % NB: Only rolling conditions  the Cartesian directions is allowed for now.
    % If duplicated  values exist, remove them.
    bc = el_bc.disp_bc;
    [bcnodes, j, i] = unique(bc.nodes);
    if(numel(bcnodes) ~= numel(bc.nodes))
        %  warning('boundary conditions have mulitiple definitions')
    end
    u_bc = bc.uu(j, :);
    mask = zeros(numel(bcnodes), G.griddim);
    % use the mask form or
    for k = 1:G.griddim
        mask(:, k) = accumarray(i, bc.mask(:, k), [numel(bcnodes), 1]);
    end
    mask = mask > 0;

    % Find the logical vector to remove Diriclet boundary conditions
    ndof = G.griddim * G.nodes.num;
    if(all(mask(:) == true))
        dirdofs = mcolon(G.griddim * (bcnodes - 1) + 1, G.griddim * (bcnodes - 1) + G.griddim)';
        u_bc = reshape(u_bc', [], 1);
    else
        dirdofs = mcolon(G.griddim * (bcnodes - 1) + 1, G.griddim * (bcnodes - 1) + G.griddim)';
        dirdofs = reshape(dirdofs, G.griddim, []);
        mm = mask';
        ind = find(mm);
        u_bc = reshape(u_bc', [], 1);
        dirdofs = reshape(dirdofs, [], 1);
        dirdofs = dirdofs(ind);
        u_bc = u_bc(ind);
    end
    isdirdofs = false(ndof, 1);
    isdirdofs(dirdofs) = true;


    %% Calulate the boundary conditions
    V_dir = nan(ndof, 1);
    V_dir(isdirdofs) = u_bc(:);
    V_dir(~isdirdofs) = 0;
    rhso = -S * V_dir;

    %% Calculate load terms
    % There are several alternatives, which may lead to different errors in particular for thin
    % long cells, see paper
    f = calculateVolumeTerm(G, load, qc_all, qcvol, opt);
    if(~isempty(opt.pressure))
        div = VEM_div(G);
        f_pressure = div'*opt.pressure;
    else
        f_pressure = zeros(size(f));
    end
    % Add load to right hand side
    rhs = rhso + f + f_pressure;
    % Add boundary forces
    if(~isempty(el_bc.force_bc))
        if(numel(el_bc.force_bc.faces) > 0)
            bc       = el_bc;
            bc_force = bc.force_bc;
            bcf      = bc_force;

            % Find weights
            faces    = bcf.faces;
            assert(numel(faces) == numel(unique(faces)));
            lnn      = G.faces.nodePos(faces + 1) - G.faces.nodePos(faces);
            inodes   = mcolon(G.faces.nodePos(bcf.faces), G.faces.nodePos(bcf.faces + 1) - 1);
            nodes    = G.faces.nodes(inodes);
            lfacenum = rldecode([1:numel(faces)]', lnn);
            assert(all(sum(G.faces.neighbors(faces, :) == 0, 2) == 1));

            if(G.griddim == 2)
                qf = qf_all(rldecode(faces(:), lnn), :);
            else
                qf = qf_all(inodes);
            end

            % Find nodes corresponding to forces
            dofs  = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim);
            force = bsxfun(@times, bc_force.force(lfacenum, :), qf);
            % Transform to degree of freedom indexing
            force = force';
            fbc   = accumarray(dofs', force(:)', [ndof, 1]);
            rhs   = rhs + fbc;
        end
    end
    % Reduce the degrees of freedom
    rhs = rhs(~isdirdofs);
    A   = S(~isdirdofs, ~isdirdofs);
    A   = (A + A') / 2; % The matrix is theoretically symmetric, make sure that it is also symmetric
                        % numerically
    x   = opt.linsolve(A, rhs);
    u   = nan(ndof, 1);

    u(isdirdofs)  = V_dir(isdirdofs);
    u(~isdirdofs) = x;
    uu = reshape(u, G.griddim, [])';

    if(nargout == 2)
        extra.A    = A;
        extra.S    = S;
        vdiv       = VEM_div(G);
        extra.disc = struct('A'         , A                        , ...
                            'isdirdofs' , isdirdofs                , ...
                            'rhs'       , rhs                      , ...
                            'V_dir'     , V_dir                    , ...
                            'ovol_div'  , vdiv                     , ...
                            'gradP'     , vdiv(:    , ~isdirdofs)' , ...
                            'div'       , vdiv(:    , ~isdirdofs)  , ...
                            'divrhs'    , vdiv * V_dir);

    end

end

function f = calculateVolumeTerm(G, load, qc_all, qcvol, opt)

    cells  = [1:G.cells.num]';
    inodes = mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells + 1) - 1');
    nodes  = G.cells.nodes(inodes);

    switch opt.force_method

      case 'node_force'
        % Evaluate forces at nodes. The result is weighted, using adjacent
        % cell volume contributions, see Paulino et al's paper (doi:10.1016/j.cma.2014.05.005)

        X = G.nodes.coords(nodes, :);
        w = qcvol;
        ll = bsxfun(@times, load(X), w)';

      case  'cell_force_baric'
        %
        % Evaluate the forces at nodes. Then, from the nodal values, compute the (exact)
        % L^2 projection on each cell and, using this cell values, assemble
        % the load term in term of the degrees of freedom.
        %
        % The VEM theory tells us that there exist virtual basis such that the
        % two steps above can be done exactly. See Ahmad et al (doi:10.1016/j.camwa.2013.05.015)
        %

        nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
        X = rldecode(G.cells.centroids(cells, :), nlc);
        lcellnum = rldecode(cells, nlc);
        BB = nan(numel(cells), 1);
        for i = 1:G.griddim
            BB(:, i) = accumarray(lcellnum, G.nodes.coords(nodes, i), [numel(cells), 1]);
        end
        fac  = accumarray(lcellnum, 1, [numel(cells), 1]);
        BB   = bsxfun(@rdivide, BB, fac);
        XB   = X - rldecode(BB, nlc);
        vols = rldecode(G.cells.volumes(cells, :), nlc);

        % weights to calculate volume force term from the nodal values.
        if(G.griddim == 3)
            w = (vols ./ rldecode(nlc, nlc) + sum(qc_all(inodes, :) .* (XB), 2));
        else
            assert(G.griddim == 2)
            w = (vols ./ rldecode(nlc, nlc) + sum(qc_all(inodes, :) .* (XB), 2));
        end
        ll   = bsxfun(@times, load(X), w)';

      case 'cell_force'
        % Evaluate the force at cell centroids. Then, for each node, sum up each
        % adjacent cell contributions after weighting them with volume
        % contribution.

        nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
        X   = rldecode(G.cells.centroids(cells, :), nlc);
        w   = qcvol;
        ll  = bsxfun(@times, load(X), w)';

      case 'dual_grad_type'
        % For the virtual basis, there exists a natural divergence operator (from node
        % values to cell values) that can be computed exactly. This operator
        % does not depend on the stabilisation term which is chosen for the
        % stiffness matrix
        %
        % By duality we can define a gradient operator (from cell-values to
        % node values). When the force can be expressed as a gradient, this
        % gives us a way to compute the load term which is implemented here.
        %
        % Such computation has appeared to be more stable, see reference.
        %
        %
        nlc     = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
        X       = rldecode(G.cells.centroids(cells, :), nlc);
        rel_vec = -(X-G.nodes.coords(nodes, :));
        ll      = bsxfun(@times, load(X), qc_all.*rel_vec)';

      otherwise
        error('No such force  calculation')
    end

    ndof = G.griddim * G.nodes.num;
    dofs = mcolon(G.griddim * (nodes - 1) + 1, G.griddim * (nodes - 1) + G.griddim)';
    f    = accumarray(dofs(:), ll(:), [ndof, 1]);

end