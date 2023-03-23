function output = runConvSim(G, params, varargin)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


    opt = struct('verbose'   , false, ...
                 'blocksize' , []   , ...
                 'useVirtual', true , ...
                 'bcetazero' , false);

    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;

    force_fun = params.force_fun;
    mu_fun = params.mu_fun;
    u_fun = params.u_fun;
    eta = params.eta;
    alpha = params.alpha;

    doVem = false;
    [tbls, mappings] = setupStandardTables(G, 'useVirtual', useVirtual);
    coltbl = tbls.coltbl;
    nodefacetbl = tbls.nodefacetbl;
    nodefacecoltbl = tbls.nodefacecoltbl;

    Nc = G.cells.num; 
    Nf = G.faces.num; 
    Nn = G.nodes.num; 
    Nd = G.griddim; 

    % prepare input for analytical functions
    for idim = 1 : Nd
        cc{idim} = G.cells.centroids(:, idim);
    end
    mu = mu_fun(cc{:});
    lambda = alpha*mu;

    prop.mu = mu; 
    prop.lambda = lambda; 

    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);

    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    bcnodefacetbl = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    bcnodefacecoltbl = crossIndexArray(bcnodefacetbl, coltbl, {}, 'optpureproduct', ...
                                       true);
    clear bcfacetbl

    [~, nodefacecents] = computeNodeFaceCentroids(G, eta, tbls, 'bcetazero', opt.bcetazero);

    map = TensorMap();
    map.fromTbl = nodefacecoltbl;
    map.toTbl = bcnodefacecoltbl;
    map.mergefds = {'nodes', 'faces', 'coldim'};
    map = map.setup();

    bcnodefacecents = map.eval(nodefacecents);
    % Here, we assume a given structure of bcnodefacecoltbl:
    bcnodefacecents = reshape(bcnodefacecents, Nd, [])';
    bcnum = bcnodefacetbl.num;

    dotest = false;
    if dotest
        % plot continuity points 
        figure
        hold on
        plotGrid(G)
        nodefacecents = reshape(nodefacecents, Nd, [])';
        if Nd == 2
            plot(nodefacecents(:, 1), nodefacecents(:, 2), '*');
        else
            plot3(nodefacecents(:, 1), nodefacecents(:, 2), nodefacecents(:, 3), '*');
        end
    end

    % Prepare input for analytical functions
    for idim = 1 : Nd
        bnfc{idim} = bcnodefacecents(:, idim);
    end

    % Compute boundary conditions
    for idim = 1 : Nd
        linform = zeros(bcnum, Nd);
        linform(:, idim) = 1;
        linforms{idim} = linform;
        linformvals{idim} = u_fun{idim}(bnfc{:});
    end

    bcfaces = bcnodefacetbl.get('faces');
    bcnodes = bcnodefacetbl.get('nodes');
    extbcnodefacetbl.faces = repmat(bcfaces, Nd, 1);
    extbcnodefacetbl.nodes = repmat(bcnodes, Nd, 1);
    extbcnodefacetbl = IndexArray(extbcnodefacetbl);

    bc.bcnodefacetbl = extbcnodefacetbl;
    bc.linform = vertcat(linforms{:});
    bc.linformvals = vertcat(linformvals{:});
    clear extbcnodefacetbl linforms linformvals

    % Compute body force
    force = NaN(Nc, Nd);
    for idim = 1 : Nd
        force(:, idim) = force_fun{idim}(cc{:});
    end
    % note minus sign (matter of convention)
    force = - bsxfun(@times, G.cells.volumes, force);
    % figure
    % quiver(cc{1}, cc{2}, force(:, 1), force(:, 2));

    % Here, we assume we know the structure of cellcoltbl;
    force = reshape(force', [], 1);

    loadstruct.bc = bc;
    loadstruct.force = force;
    loadstruct.extforce = zeros(tbls.nodefacecoltbl.num, 1);
    clear bc force

    if ~isempty(opt.blocksize)
        assembly = blockAssembleMPSA(G, prop, loadstruct, eta, tbls, mappings, ...
                                     'blocksize', opt.blocksize, 'verbose', ...
                                     true, 'useVirtual', useVirtual);
    else
        assembly = assembleMPSA(G, prop, loadstruct, eta, tbls, mappings, ...
                                'bcetazero', opt.bcetazero);
    end

    clear prop loadstruct

    B   = assembly.B  ;
    rhs = assembly.rhs;
    sol = B\rhs;

    % Displacement values at cell centers.
    cellcoltbl = tbls.cellcoltbl;
    n = cellcoltbl.num;

    u = sol(1 : n);
    u = reshape(u, Nd, [])';    

    output = struct('B'       , B       , ...
                    'assembly', assembly, ...
                    'rhs'     , rhs     , ...
                    'tbls'    , tbls    , ...
                    'u'       , u);
end
