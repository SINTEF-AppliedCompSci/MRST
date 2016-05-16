function p = partitionFracture(G, pm, nw, varargin)
% Generates partition vector for the fracture grid. See getRsbGridsHFM for
% description of input parameters.

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


if isstruct(varargin{1})
    opts = varargin{1};
else
    opts = struct( 'partition_frac'    , true   , ...
                   'dof_frac'          , 1      , ...
                   'sysMatrix'         , []     , ...
                   'use_metisF'        , true   , ...
                   'coarseDimsF'       , []     , ...
                   'sampleDimsF'       , []     , ...
                   'paddedPartitionF'  , false  );
    opts = merge_options(opts, varargin{:});
end
p = zeros(G.cells.num,1);
p(1:numel(pm)) = pm;
frac_cells = G.Matrix.cells.num+1:G.cells.num;
if opts.partition_frac == true
    if opts.use_metisF
        if isempty(opts.sysMatrix)
            fprintf('\nSystem matrix not provided, using unweighted adjacency graph to partition fractures.\n');
            A = getConnectivityMatrix(getNeighbourship(G,'Topological'));
        else
            A = opts.sysMatrix;
        end
        Af = A(frac_cells,frac_cells);
        Af = Af - diag(diag(Af));
        Af = Af + abs(diag(sum(Af,2)));
        mopt = struct('ufactor', 1.5, ...
            'ncuts',   10, ...
            'seed',    0, ...
            'useLog',  false, ...
            'no2hop',  false);
        metis_options = sprintf('-minconn -contig -ufactor=%d -objtype=cut -ncuts=%d -seed=%d', ...
            floor(10*(mopt.ufactor - 1)), mopt.ncuts, mopt.seed);
        if mopt.no2hop
            metis_options = [metis_options, ' -no2hop'];
        end
        for i = 1:numel(nw)
            nwcells = [];
            for j = 1:numel(nw(i).lines);
                ll = nw(i).lines(j);
                start = G.FracGrid.(['Frac',num2str(ll)]).cells.start;
                finish = G.FracGrid.(['Frac',num2str(ll)]).cells.start + ...
                    G.FracGrid.(['Frac',num2str(ll)]).cells.num-1;
                nwcells = [nwcells,start:finish]; %#ok
            end
            if opts.dof_frac(i) == 1
                p(nwcells) = max(p) + 1;
            elseif opts.dof_frac(i) == numel(nwcells)
                p(nwcells) = max(p) + transpose(1:numel(nwcells));
            else
                Afn = Af(nwcells-G.Matrix.cells.num,nwcells-G.Matrix.cells.num); % only considering fracture-fracture connectivity inside a network
                pf = callMetisMatrix(Afn, opts.dof_frac(i), metis_options);
                pf = compressPartition(pf);
                p(nwcells) = pf+max(p);
            end
        end
    else
        lset = {}; count = 1;
        for i = 1:numel(nw)
            for j = 1:numel(nw(i).lines);
                ll = nw(i).lines(j);
                F = G.FracGrid.(['Frac',num2str(ll)]);
                start = F.cells.start;
                finish = F.cells.start + ...
                    F.cells.num-1;
                lset{count,1} = start:finish; %#ok
                Gf = extractSubgrid(G,lset{count,1});
                Gf.type = [F.type,Gf.type];
                Gf.layerSize = F.layerSize;
                Gf.numLayers = F.numLayers;
                Gf.griddim = Gf.griddim - 1;
                pf = getPartitionVectorFracPlane(Gf,F,opts);
                p(lset{count,1}) = pf+max(p);
                count = count + 1;
            end
        end
    end
else
    % no coarsening/multiscale in fractures. Fracture fine-cells belong to
    % matrix coarse cells
    warning(sprintf(['''partition_frac'' key is set to false in function arguments.',...
        '\nNo multiscale formulation in fractures.']));%#ok
    for i = 1:numel(frac_cells)
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        matrix_cells = nnc(nnc(:,2)==frac_cells(i),1);
        coarse_p = unique(pm(matrix_cells));
        p(frac_cells(i)) = coarse_p(1);
    end
end
return


function pf = getPartitionVectorFracPlane(Gf,F,opts)
dims = max(Gf.nodes.coords);
if isempty(opts.sampleDimsF)
    if strcmp(Gf.type{1,1},'tensorGrid') || strcmp(Gf.type{1,1},'cartGrid')
        if isfield(F,'cartDims')
            Gf.cartDims = F.cartDims;
        elseif any(strcmp(Gf.type,'layered'))
            Gf.cartDims = [F.layerSize F.numLayers];
        end
        if opts.paddedPartitionF
            pf = partitionUniformPadded(Gf, opts.coarseDimsF);
        else
            pf = partitionUI(Gf, opts.coarseDimsF);
        end
    end
else
    assert(~isempty(opts.sampleDimsF),'Provide sampling grid dimensions using ''sampleDims'' key');
    Gf = rmfield(Gf,'cartDims');
    Gc = cartGrid(opts.sampleDimsF, dims);
    % Make coarse grid
    if opts.paddedPartitionF
        pf = partitionUniformPadded(Gc, opts.coarseDimsF);
    else
        pf = partitionUI(Gc, opts.coarseDimsF);
    end
    pf = reshape(pf, Gc.cartDims);
    Gf.griddim = Gf.griddim + 1;
    pf = sampleFromBox(Gf, pf);
end
pf = processPartition(Gf,compressPartition(pf));
return