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
                 'useVirtual', false, ...
                 'bcetazero' , false);

    opt = merge_options(opt, varargin{:});

    useVirtual = opt.useVirtual;

    force_fun = params.force_fun;
    mu_fun    = params.mu_fun;
    u_fun     = params.u_fun;
    eta       = params.eta;
    alpha     = params.alpha;

    [tbls, mappings] = setupMpsaStandardTables(G, 'useVirtual', useVirtual);
    vectbl         = tbls.vectbl;
    nodefacetbl    = tbls.nodefacetbl;
    nodefacevectbl = tbls.nodefacevectbl;

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

    prop.mu     = mu; 
    prop.lambda = lambda; 

    isBoundary = any(G.faces.neighbors == 0, 2); 
    bcfaces =  find(isBoundary);

    bcfacetbl.faces = bcfaces;
    bcfacetbl = IndexArray(bcfacetbl);
    [bcnodefacetbl, indstruct] = crossIndexArray(bcfacetbl, nodefacetbl, {'faces'});
    if useVirtual
        bcface_from_bcnodeface   = indstruct{1}.inds;
        nodeface_from_bcnodeface = indstruct{2}.inds;
    else
    end
    bcnodefacevectbl = crossIndexArray(bcnodefacetbl, vectbl, {}, ...
                                       'optpureproduct', true   , ...
                                       'virtual', useVirtual);
    clear bcfacetbl

    [~, nodefacecents] = computeNodeFaceCentroids2(G, eta, tbls, mappings,  ...
                                                   'bcetazero', opt.bcetazero, ...
                                                   'useVirtual', useVirtual);

    map = TensorMap();
    map.fromTbl  = nodefacevectbl;
    map.toTbl    = bcnodefacevectbl;
    map.mergefds = {'nodes', 'faces', 'vec'};

    if useVirtual
        map.pivottbl = bcnodefacevectbl;

        [vec, i] = ind2sub([vectbl.num, bcnodefacetbl.num], (1 : bcnodefacevectbl.num)');
        map.dispind1 = sub2ind([vectbl.num, nodefacetbl.num], vec, nodeface_from_bcnodeface(i));
        map.dispind2 = (1 : bcnodefacevectbl.num)';
        map.issetup = true;
        
    else
        map = map.setup();
    end
        
    bcnodefacecents = map.eval(nodefacecents);
    % Here, we assume a given structure of bcnodefacevectbl:
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

    % Here, we assume we know the structure of cellvectbl;
    force = reshape(force', [], 1);

    loadstruct.bc = bc;
    loadstruct.force = force;
    loadstruct.extforce = zeros(tbls.nodefacevectbl.num, 1);
    clear bc force

    if ~isempty(opt.blocksize)
        assembly = blockAssembleMPSA2(G, prop, loadstruct, eta, tbls, mappings, ...
                                      'blocksize' , opt.blocksize, ...
                                      'verbose'   , true         , ...
                                      'useVirtual', useVirtual);
    else
        assembly = assembleMPSA2(G, prop, loadstruct, eta, tbls, mappings, ...
                                 'bcetazero', opt.bcetazero, ...
                                 'useVirtual', useVirtual);
    end

    clear prop loadstruct

    B   = assembly.B  ;
    rhs = assembly.rhs;
    sol = B\rhs;

    % Displacement values at cell centers.
    cellvectbl = tbls.cellvectbl;
    n = cellvectbl.num;

    u = sol(1 : n);
    u = reshape(u, Nd, [])';    

    output = struct('B'       , B       , ...
                    'assembly', assembly, ...
                    'rhs'     , rhs     , ...
                    'tbls'    , tbls    , ...
                    'u'       , u);
end
