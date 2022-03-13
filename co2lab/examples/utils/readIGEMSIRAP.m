function G = readIGEMSIRAP(name, i, varargin)
%Read and process IGEMS .irap grids.
%
% SYNOPSIS:
%   G = readIGEMSIRAP('flatUP1', 3)
%   G = readIGEMSIRAP('flatUP1', 3, 'coarse', [5 5])
%   G = readIGEMSIRAP('1_FMM.irap', [], 'dir', idir, 'save', false);
%
% PARAMETERS:
%   name   - The name of the IGEMS grid. For valid names, run
%            >>> ls(fullfile(getDatasetPath('IGEMSsurfaces'), 'surfaces'))
%
%   i      - The index of the realization of the grid. A value from 1 to
%            100.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              coarse -- Vector of 2 elements corresponding to a coarsening
%              factor for logical i and j directions respectively. A coarse
%              vector of [2, 2] will produce a grid 1/2^2 = 1/4th the size
%              of [1, 1] and so on.
%
%              save   -- Bolean. If true, a mat file with the constructed
%              grid is saved for later use. Default: false
%
%              dir    -- String. Load files from this directory.
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
%
% SEE ALSO:
%   `readIrapClassicAsciiSurf`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('coarse', [1 1], ...
             'save', false, ...
             'dir',  []);
opt = merge_options(opt, varargin{:});
coarse = opt.coarse;

subset = false;
if isempty(opt.dir)
   %result_dir = fullfile(mrstPath('co2lab'), 'data', 'igems');
   info = getDatasetInfo('IGEMSsurfaces');
   if i == 1 && ~exist(fullfile(mrstDataDirectory(), info.name), 'dir')
       file_dir = fullfile(getDatasetPath('IGEMSsample'), 'one_of_each');
       subset = true;
   else
       file_dir = fullfile(getDatasetPath('IGEMSsurfaces'), 'surfaces');
   end
else
   result_dir = opt.dir;
   file_dir   = opt.dir;
end

if all(coarse == 1)
    outfilename = ['grid', name];
else
    outfilename = ['grid', name, '_', sprintf('%d', coarse)];
end

if ~isempty(i)
    outfilename = [outfilename, '_', num2str(i)];
end

outfilename = [outfilename, '.mat'];

try
    load(outfilename)
catch %#ok<*CTCH>
    if isempty(strfind(name,'.irap')) && isempty(i)
       name = [name '.irap'];
    end
    if(isempty(i))
       filename=[file_dir, filesep, name];
    else
        if subset
            filename=[file_dir, filesep, int2str(i),'_', name '.irap'];
        else
            filename=[file_dir, filesep, name, filesep, int2str(i),'_', name '.irap'];
        end
    end
    [x, y, Z]= readIrapClassicAsciiSurf(filename);
    assert(all(diff(x)==100));
    assert(all(diff(y)==100));
    if ~(prod(coarse)==1)
        Z=Z(1:coarse(1):end,1:coarse(2):end);
        G=cartGrid([size(Z)-1,1],[(size(Z)-1).*coarse*100,100]);
    else
        G=cartGrid([size(Z)-1,1],[(size(Z)-1),1]*100);
    end
    Z(:,1)   = Z(:,2);
    Z(:,end) = Z(:,end-1);
    nxy = prod(G.cartDims(1:2)+1);
    G.nodes.coords(1:nxy,3) = G.nodes.coords(1:nxy,3)+Z(:);
    G.nodes.coords(nxy+1:end,3) = G.nodes.coords(nxy+1:end,3)+Z(:);
    Z_cells =  (Z(1:G.cartDims(1),1:G.cartDims(2))+...
                Z(2:end,1:G.cartDims(2))+...
                Z(1:G.cartDims(1),2:end)+...
                Z(2:end,2:end))/4;
    assert(numel(Z_cells)==G.cells.num);
    
    rmcells=find(~isfinite(Z_cells(:)));
    G=removeCells(G,rmcells);
    try
        G = mcomputeGeometry(G);
    catch
        G = computeGeometry(G);
    end


    if opt.save,
       mkdirif([result_dir,filesep, name])
       mkdirif([result_dir,filesep, name, filesep, 'grid'])
       save(outfilename,'-v7.3','G')
    end
end
end

function mkdirif(directory)
   if ~exist(directory, 'dir')
       mkdir(directory);
   end
end
