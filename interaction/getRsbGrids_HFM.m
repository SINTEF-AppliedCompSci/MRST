function [CG, CGf] = getRsbGrids_HFM(G, F, nw, varargin)
% Computes coarse grid and interaction region for a grid with fractures
%
% SYNOPSIS:
%   [CG, CGf] = getRsbGrids_HFM(G, F, nw)
%   [CG, CGf] = getRsbGrids_HFM(G, F, nw, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function coarsens the matrix and fracture grid separately and
%   computes supports for the fracture and matrix basis functions. All of
%   this is then combined into one coarse grid to be used for the
%   multiscale solver. There are no restrictions on grid definition.
%   Fracture partitioning algorithms are graph based.
%
% REQUIRED PARAMETERS:
%
%   G  - Grid structure with fractures as defined by assembleGlobalGrid.
%
%   F  - Structure containing information about partitioned fracture
%        lines as returned by assembleFracNodes2D. See also gridFracture.
%
%   nw - same as fracture.networks as returned by getIndepNetwork. Contains
%        line indices for independant fracture networks
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   coarsen     - Number of fine cells per coarse block in each physical
%                 direction used to extract the desirable coarse grid
%                 resolution (See partitionUI). Assumed to be a LENGTH 2 or
%                 3 vector of integers. Used only for cartesian meshes if
%                 use_metis == false.
%
%   paddedPartition - If true, uses partitionUniformPadded to partition
%                     matrix coarse grid. Requires 'coarsen'.
%
%   use_metis   - Use metis to partition matrix.
%
%   dof_matrix  - Degrees of freedom in matrix at coarse scale. Required is
%                 use_metis == true.
%
%   MatrixTrans - Matrix transmissibility. Required if use_metis == true.
%
%   partition_frac - If true, fractures will be partitioned and there will
%                    be separate interaction regions for fractures.
%                    Otherwise fractures will be a part of the matrix
%                    coarse grid and interaction regions.
%
%   dof_frac    - Degrees of freedom in each fracture network at coarse
%                 scale. Must either be a single value < minimum number of
%                 fine cells in all independant fracture networks or 1
%                 value per intependant fracture network.
%
%   sysMatrix   - Fine scale system (transmissibility) matrix A (see
%                 incompTPFA). Must be specified if use_metis or
%                 partition_frac are true
%
%   fracIRtype  - Matrix support type for the fracture basis. See
%                 storeInteractionRegionFrac
%
% RETURNS:
%   CG - Coarse grid for grid G with basis supports and coarse node
%        indices. See storeInteractionRegion, generateCoarseGrid.
%
%   CGf - Fracture coarse grid with similar information as CG but only for
%         the internal fracture grid Gf as returned by assembleFracGrid.
%
% SEE ALSO:
%   partitionUI, partitionUniformPadded, storeInteractionRegion,
%   storeInteractionRegionFrac, partitionMETIS, callMetisMatrix.

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


opts = struct(...
    'coarsen'         , []     , ...
    'paddedPartition' , true   , ...
    'use_metis'       , false  , ...
    'dof_matrix'      , 5      , ...
    'MatrixTrans'     , computeTrans(G.Matrix,G.Matrix.rock), ...
    'partition_frac'  , true   , ...
    'dof_frac'        , 1      , ...
    'sysMatrix'       , []     , ...
    'fracIRtype'      , 4      );

opts = merge_options(opts, varargin{:});
Gm = G.Matrix; nt = G.cells.num;

if opts.partition_frac == true
    assert(all(opts.dof_frac>=1),['Degrees of freedom for each fracture ',...
        'network at coarse scale must be passed as an integer >= 1.']);
    assert(~isempty(opts.sysMatrix) && isequal(size(opts.sysMatrix),[nt nt]),...
        sprintf(['System matrix for the given grid ''G'' is not specified.\n',...
        'It is needed for partitioning fractures using METIS.']));
    numfc = zeros(1,numel(F));
    numnc = zeros(1,numel(nw));
    for i = 1:numel(nw)
        for j = 1:numel(nw(i).lines)
            numfc(nw(i).lines(j)) = F(nw(i).lines(j)).cells.num;
        end
        numnc(i) = sum(numfc(nw(i).lines));
    end
    assert(numel(opts.dof_frac)==1 || numel(opts.dof_frac)==numel(nw),...
        sprintf(['Specify either 1 DOF per fracture network or 1 DOF for ',...
        'each fracture network.\nThere are ',num2str(numel(nw)),' independant fracture networks.']));
    if numel(opts.dof_frac)==1 && numel(opts.dof_frac)<numel(nw) 
        warning(sprintf(['It would be a good idea to specify 1 coarse DOF ',...
            ' value per fracture network.\nProceeding with the a coarse ',...
            ' DOF value of ',num2str(opts.dof_frac),' for each fracture network.'])); %#ok
        [maxdof,locmin] = min(numfc);
        assert(opts.dof_frac<=maxdof,sprintf(['Coarse degrees of freedom cannot be ',...
            'more than fine-scale degrees of freedom.\nLine ',num2str(locmin),...
            ' has ',num2str(maxdof),' fine cells and coarse DOF specified is ',...
            num2str(opts.dof_frac),' per fracture line.']));     
        opts.dof_frac = repmat(opts.dof_frac,numel(nw),1);
    else
        assert(all(numnc>=opts.dof_frac),sprintf(['Coarse DOF per fracture ',...
            'network cannot be more than fine-scale DOF per fracture network.\n',...
            'Reduce entries ',num2str(transpose(find(numnc<opts.dof_frac))),...
            ' of ''dof_frac'' argument.']));
    end
end
if opts.use_metis==true
    assert(opts.dof_matrix>1,sprintf(['Number of blocks into which to ',...
        'partition the grid ''G'' is not specified.\nSpecify a ',...
        'positive integer in the function call using the ''dof_matrix'' flag.']));
    assert(~isempty(opts.MatrixTrans),sprintf(['Forgot to pass transmissibility ',...
        'vector for matrix cells ?\nAdd using key ''MatrixTrans'' in argument list.']));
else
    assert(isfield(Gm,'cartDims'),'Underlying matrix grid must be cartesian if partitionUI is to be used');
end
if opts.use_metis
    % partition only matrix grid first
    % pass metis options as 'key'/value pairs. 'useLog' may be useful here
    dof_matrix = opts.dof_matrix;
    pm = partitionMETIS(Gm, opts.MatrixTrans, dof_matrix); 
else
    coarsedims = ceil(Gm.cartDims./opts.coarsen);
    % Generate partition vector
    if opts.paddedPartition
        pm = partitionUniformPadded(Gm, coarsedims);
    else
        pm = partitionUI(Gm, coarsedims);
    end
end
p = partitionFrac(G, pm, F, nw, opts); % Assign p to fractures or partition them separately
pf = p(G.Matrix.cells.num+1:end)-max(p(1:G.Matrix.cells.num));

% Coarse grid structure from partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);

if opts.partition_frac==true
    Gf = assembleFracGrid(G);
    CGm = generateCoarseGrid(Gm,pm);
    CGf = generateCoarseGrid(Gf, pf);
    CGm = coarsenGeometry(CGm);
    CGf = coarsenGeometry(CGf);
    CGm = storeInteractionRegion(CGm);
    [CG,CGf] = storeInteractionRegionFrac(CG,CGf,CGm,...
        'type',opts.fracIRtype,'sysMatrix',opts.sysMatrix);
else
    CG = storeInteractionRegion(CG);
end
return

%-------------------------------------------------------------------------%

function pnew = partitionFrac(G, p, F, nw, opts)
pnew = zeros(G.cells.num,1);
pnew(1:numel(p)) = p;
frac_cells = G.Matrix.cells.num+1:G.cells.num;
if opts.partition_frac == true
    Af = opts.sysMatrix(G.Matrix.cells.num+1:end,G.Matrix.cells.num+1:end);
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
            nwcells = [nwcells,F(ll).cells.start:(F(ll).cells.start+F(ll).cells.num-1)]; %#ok
        end
        if opts.dof_frac(i) == 1
            pnew(nwcells) = max(pnew) + 1;
        else
            Afn = Af(nwcells-G.Matrix.cells.num,nwcells-G.Matrix.cells.num); % only considering fracture-fracture connectivity inside a network
            pf = callMetisMatrix(Afn, opts.dof_frac(i), metis_options);
            pf = compressPartition(pf);
            pnew(nwcells) = pf+max(pnew);
        end
    end
else
    % no coarsening/multiscale in fractures. Fracture fine-cells belong to
    % matrix coarse cells
    warning(sprintf(['''partition_frac key'' is set to false in function arguments.',...
        '\nNo multiscale formulation in fractures i.e. no rsb in fractures.']));%#ok
    for i = 1:numel(frac_cells)
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        matrix_cells = nnc(nnc(:,2)==frac_cells(i),1);
        coarse_p = unique(p(matrix_cells));
        pnew(frac_cells(i)) = coarse_p(1);
    end
end
return