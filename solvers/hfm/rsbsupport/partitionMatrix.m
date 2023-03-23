function pm = partitionMatrix(G,varargin)
% Generates partition vector for the matrix grid. See getRsbGridsHFM for
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
    opts = struct(  'coarseDims'                      , []     , ...
                    'sampleDims'                      , []     , ...
                    'paddedPartition'                 , false  , ...
                    'use_metis'                       , false  , ...
                    'dof_matrix'                      , 5      , ...
                    'matrixTrans'                     , computeTrans(G.Matrix,G.Matrix.rock), ...
                    'Wells'                           , []     , ...
                    'nearWellRefinement'              , false  , ...
                    'nearWellRefinementMultipler'     , 5      );
    opts = merge_options(opts, varargin{:});
end
Gm = G.Matrix;
if opts.use_metis == false && ~(isfield(Gm,'cartDims'))
    fprintf('\n Underlying matrix grid is not cartesian. Sampling partition from box.\n');
end

if opts.use_metis
    pm = getPartitionVectorMETIS(Gm,opts);
else
    pm = getPartitionVectorUI(G,Gm,opts);
end
return

%-------------------------------------------------------------------------%

function p = getPartitionVectorMETIS(Gm,opts)
dof = opts.dof_matrix;
if ~opts.nearWellRefinement
    p = partitionMETIS(Gm, opts.matrixTrans, dof);
    return
else
    T = opts.matrixTrans;
    multiplier = opts.nearWellRefinementMultipler;
    wcells = [opts.Wells.cells];
    dof_i = dof*multiplier*(dof*multiplier<=ceil(Gm.cells.num/multiplier)) + ...
        mean(Gm.cells.num,dof)*(dof*multiplier>ceil(Gm.cells.num/multiplier));
    p0 = transpose(1:Gm.cells.num);
    p1 = partitionMETIS(Gm,T,dof_i);
    pm = partitionMETIS(Gm,T,dof);
    p = pm;
    for i = 1:numel(wcells)
        wcell = wcells(i);
        
        sub1 = pm == pm(wcell); % Fine cells of primary coarse cell containing well
        p(sub1) = max(p) + p1(sub1);
        
        sub2 = p1 == p1(wcell);
        p(sub2) = max(p) + p0(sub2);
    end
end
p = processPartition(Gm,compressPartition(p));
return

%-------------------------------------------------------------------------%

function p = getPartitionVectorUI(G,Gm,opts)
if ~opts.nearWellRefinement
    p = getPartitionVector(G,Gm,opts);
    return
else
    multiplier = opts.nearWellRefinementMultipler;
    wcells = [opts.Wells.cells];
    opts2 = opts; % opts3 = opts;
    opts2.coarseDims = floor(opts.coarseDims.*multiplier);
%     opts3.coarseDims = floor(opts.coarseDims.*multiplier*4);
    p0 = transpose(1:Gm.cells.num);
    p1 = getPartitionVector(G,Gm,opts2);
%     p2 = getPartitionVector(G,Gm,opts3);
    pm = getPartitionVector(G,Gm,opts);
    p = pm;
    for i = 1:numel(wcells)
        wcell = wcells(i);
        sub1 = pm == pm(wcell); % Fine cells of primary coarse cell containing well
        p(sub1) = max(p) + p1(sub1);
        
        sub2 = p1 == p1(wcell);
        p(sub2) = max(p) + p0(sub2);
        
%         sub3 = p2 == p2(wcell);
%         p(sub3) = max(p) + p0(sub3);
    end
end
p = processPartition(Gm,compressPartition(p));
return

%-------------------------------------------------------------------------%

function pm = getPartitionVector(G,Gm,opts)
dims = max(G.nodes.coords);
if isfield(Gm,'cartDims')
    % Generate partition vector
    if strcmp(G.type{1,1},'processGRDECL')
        if isempty(opts.sampleDims), opts.sampleDims = G.cartDims; end
        Gc = cartGrid(opts.sampleDims, dims);
        % Make coarse grid
        pm = partitionUI(Gc, opts.coarseDims);
        pm = reshape(pm, Gc.cartDims);
        pm = sampleFromBox(Gm, pm);
    elseif opts.paddedPartition
        pm = partitionUniformPadded(Gm, opts.coarseDims);
    else
        pm = partitionUI(Gm, opts.coarseDims);
    end
else
    assert(~isempty(opts.sampleDims),'Provide sampling grid dimensions using ''sampleDims'' key');
    Gc = cartGrid(opts.sampleDims, dims);
    % Make coarse grid
    if opts.paddedPartition
        pm = partitionUniformPadded(Gc, opts.coarseDims);
    else
        pm = partitionUI(Gc, opts.coarseDims);
    end
    pm = reshape(pm, Gc.cartDims);
    pm = sampleFromBox(Gm, pm);
end
pm = processPartition(Gm,compressPartition(pm));
return