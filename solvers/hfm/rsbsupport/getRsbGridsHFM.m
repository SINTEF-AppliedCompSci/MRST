function [CG, CGf] = getRsbGridsHFM(G, nw, varargin)
% Computes coarse grid and interaction region for a grid with fractures
%
% SYNOPSIS:
%   [CG, CGf] = getRsbGridsHFM(G, nw)
%   [CG, CGf] = getRsbGridsHFM(G, nw, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function coarsens the matrix and fracture grids separately and
%   computes support regions for the fracture and matrix basis functions.
%   All of this is then combined into one coarse grid to be used for the
%   multiscale solver. There are no restrictions on grid definition.
%   Fracture partitioning algorithms are graph based.
%
% REQUIRED PARAMETERS:
%
%   G  - Grid structure with fractures as defined by assembleGlobalGrid.
%
%   nw - Structure containing indices of fractures that make up every
%        independent fracture network. same as fracture.networks for 2D
%        examples (see getIndepNetwork)
%
% OPTIONAL PARAMETERS:
%
%   pm - Partition vector mapping matrix fine cells to coarse blocks.
%        Size = nm-by-1 where nm is the total number of fine cells in the
%        matrix.
%
%   coarseDims  -  Number of coarse blocks in each physical direction.
%                  Assumed to be a LENGTH 2 or 3 vector of integers. To be
%                  used only when the matrix grid is structured.
%
%   sampleDims  - Hypothetical cartesian dimensions (specifying number of
%                 fine-cells in each coordinate direction) representing a
%                 unstructured grid. Used to construct a structured coarse
%                 grid, in the same physical space, which is then utilized
%                 as a sample for mapping the unstructured fine grid to the
%                 structured coarse grid. See sampleFromBox.
%
%   use_metis   - Logical value (true or false) to indicate the use of
%                 metis to partition matrix.
%
%   dof_matrix  - Degrees of freedom in matrix at coarse scale. Required if
%                 use_metis = true.
%
%   MatrixTrans - Transmissibility vector for all fine-grid faces in the
%                 matrix
%
%   partition_frac - Logical value to indicate if fractures must be
%                    coarsened. Otherwise fractures will be a part of the
%                    matrix coarse grid and support regions and there will
%                    be no coarse nodes inside a fracture. Useful to
%                    observe the impact of coarsening fractures.
%
%   coarseDimsF - Number of coarse blocks for a fracture plate in each
%                   physical direction. Aplicable only to rectangular
%                   fracture plates in 3D system.
%
%   sampleDimsF - Similar to its counterpart in the matrix. Aplicable only
%                 to rectangular fracture plates in 3D system.
%
%   use_metisF  - Logical value (true or false) to indicate the use of
%                 metis to partition fractures.
%
%   dof_frac     - Degrees of freedom in each fracture network at coarse
%                  scale. Required when use_metisF = true. Must either be a
%                  single value < minimum number of fine cells in all
%                  independant fracture networks or a vector containing a
%                  value for each intependent fracture network.
%
%   sysMatrix    - Fine scale system (transmissibility) matrix A (see
%                  incompTPFA). Not necessary, but if available, allows the
%                  use of an existing variable to construct an Adjacency
%                  (or connectivity) matrix.
%
%   fracSupportRadius - levels of connectivity used to define the support
%                       region for a fracture coarse block. In other words,
%                       the topological radius of the fracture support
%                       region. See storeFractureInteractionRegion.
%
%   fullyCoupled - Coupled fracture and matrix basis functions. Allows
%                  matrix support to extend into fractures. Increases
%                  accuracy but reduces speed.
%
%   excludeBoundary - Excludes fine cells on the boundary from support
%                     regions of internal coarse blocks. Improves accuracy,
%                     especially when boundary conditions are specified.
%
%
%%%%%%%%%%%%%%%%%%%-----------EXPERIMENTAL-----------%%%%%%%%%%%%%%%%%%%%%%
%
%
%   removeCenters    - Logical value to essentially enforce a basis
%                      function value of 1 at each coarse node. May improve
%                      convergence rate in an iterative multiscale
%                      strategy.
%
%   coarseNodeOption - Possible methods for optimizing coarse node location
%                      inside fractures.
%
%   paddedPartition  - If true, uses partitionUniformPadded to partition
%                      matrix coarse grid. Requires 'coarseDims'.
%
%   paddedPartitionF - If true, uses partitionUniformPadded to partition
%                      fracture coarse grid. Requires 'coarseDimsF'.
% 
%   Wells            - Well structure. See addWell.
%
%   nearWellRefinement - Refines the coarse grid near wells by recursively
%                        splitting the coarse block containing a well a
%                        specified number of times.
%
%   nearWellRefinementMultipler - Integer denoting the number of times the
%                                 coarse block containing a well must be
%                                 recursively split into 2.
%
%
% RETURNS:
%   CG  - Coarse grid for grid G with basis supports and coarse node
%         indices. See storeInteractionRegion, generateCoarseGrid.
%
%   CGf - Fracture coarse grid with similar information as CG but only for
%         the internal fracture grid Gf as returned by assembleFracGrid.
%
% SEE ALSO:
%   partitionMatrix, partitionFracture, generateCoarseGrid,
%   storeInteractionRegion, storeInteractionRegionFrac, partitionMETIS,
%   callMetisMatrix.

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


opts = struct('pm'                       , []     , ... 
              'coarseDims'               , []     , ...
              'sampleDims'               , []     , ...
              'use_metis'                , false  , ...
              'dof_matrix'               , 5      , ...
              'matrixTrans'              , computeTrans(G.Matrix,G.Matrix.rock), ...
              'partition_frac'           , true   , ...
              'use_metisF'               , true  , ...
              'coarseDimsF'              , []     , ...
              'sampleDimsF'              , []     , ...
              'dof_frac'                 , 1      , ...
              'sysMatrix'                , []     , ...
              'fracSupportRadius'        , 8      , ...
              'fullyCoupled'             , true   , ...
              'excludeBoundary'          , true   , ...
              'removeCenters'            , true   , ...
              'coarseNodeOption'         , 'useFineCellCentroids', ...
              'paddedPartition'          , false  , ...
              'paddedPartitionF'         , false  , ...
              'Wells'                    , []     , ...
              'nearWellRefinement'       , false  , ...
              'nearWellRefinementMultipler'  , 2  );

opts = merge_options(opts, varargin{:});
opts = check_input_varargin(G,nw,opts); pm = opts.pm;
%-------------------------------------------------------------------------%
if isempty(pm), pm = partitionMatrix(G, opts); end
p  = partitionFracture(G, pm, nw, opts); % Assign p to fractures or partition them separately
% p = processPartition(G,compressPartition(p));
pf = p(G.Matrix.cells.num+1:end)-max(p(1:G.Matrix.cells.num));

% Coarse grid structure from partition vector
CG = generateCoarseGrid(G, p);
% Add centroids / geometry information on coarse grid
CG = coarsenGeometry(CG);

if opts.partition_frac==true
    CGm = getRsbGridsMatrix(G,pm,opts);
    Gf = assembleFracGrid(G);
    CGf = generateCoarseGrid(Gf, pf);
    CGf = coarsenGeometry(CGf);
    [CG,CGf] = storeFractureInteractionRegion(CG, CGf, CGm, ...
               'fullyCoupled', opts.fullyCoupled, ...
               'excludeBoundary', opts.excludeBoundary, ...
               'removeCenters', opts.removeCenters, ...
               'coarseNodeOption', opts.coarseNodeOption, ...
               'levels', opts.fracSupportRadius);
else
    CGf = struct;
    CG = storeInteractionRegion(CG);
end
return

%-------------------------------------------------------------------------%

function opts = check_input_varargin(G,nw,opts)
if opts.partition_frac == true
    numfc = zeros(1,numel(fieldnames(G.FracGrid)));
    numnc = zeros(1,numel(nw));
    for i = 1:numel(nw)
        for j = 1:numel(nw(i).lines)
            numfc(nw(i).lines(j)) = G.FracGrid.(['Frac',num2str(nw(i).lines(j))]).cells.num;
        end
        numnc(i) = sum(numfc(nw(i).lines));
    end
    % Frac DOF Related Warnings
    assert(all(opts.dof_frac>=1),['Degrees of freedom for each fracture ',...
        'network at coarse scale must be passed as an integer >= 1.']);
    assert(numel(opts.dof_frac)==1 || numel(opts.dof_frac)==numel(nw),...
        sprintf(['Specify either 1 DOF per fracture network or 1 DOF for ',...
        'each fracture network.\nThere are ',num2str(numel(nw)),' independant fracture networks.']));
    if numel(opts.dof_frac)==1 && numel(opts.dof_frac)<numel(nw)
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

if opts.nearWellRefinement
    assert(~isempty(opts.Wells),'Provide cells penetrated by wells');
end
return