function G = readIGEMSIRAP(name, i, varargin)
%Read and process IGEMS .irap grids.
%
% SYNOPSIS:
%   G = readIGEMSIRAP('flatUP1', 3)
%   G = readIGEMSIRAP('flatUP1', 3, 'coarse', [5 5])
%
% PARAMETERS:
%   name   - The name of the IGEMS grid. For valid names, run
%            >>> ls(fullfile(VEROOTDIR, 'data', 'igems', 'surfaces'))
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
%              grid is saved for later use. Default: true
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
%
% SEE ALSO:
%   readIrapClassicAsciiSurf

%{
#COPYRIGHT#
%}

require mex
opt = struct('coarse', [1 1], 'save', true);
opt = merge_options(opt, varargin{:});
coarse = opt.coarse;


result_dir = fullfile(VEROOTDIR, 'data', 'igems');
file_dir = fullfile(result_dir, 'surfaces');

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
catch
    if(isempty(i))
       filename=[file_dir, filesep, name,'.irap'];
    else
       filename=[file_dir, filesep ,name, filesep, int2str(i),'_',name,'.irap'];
    end
    [x, y, Z, angle]= readIrapClassicAsciiSurf(filename);
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
