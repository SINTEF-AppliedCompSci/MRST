classdef CellDataPatch < handle
    properties
        fastMode = false % is set to true, only slices are plotted, not lids
        logScale  = false;
        patchOpt  = {};
    end
    properties (SetAccess = protected)
        patchMain      % standard patch (for non-play mode)
        patchSlices    % playmode - patches
        patchLids      % playmode - patches
        main_ix   = [] % used for resetting color/filter
        slice_ix  = [] % used for resetting color/filter
        lid_ix    = [] % used for resetting color/filter
        patchType      % 'single', 'multiple' 
        binMeans  = [] % means of bins 
        sliceThreshold = 10;
    end
    properties (SetAccess = immutable)
        nbins          % number of bins used for play-mode        
    end
    properties (Hidden)
        gridData       % store all data needed for changing cell-set/filter-data/color-data
    end
    properties (Dependent) % set/get - props
        cells          % logic index to cells
        colorData
        filterData
        value
        curBin         % current bin
        playMode       % true/false
        Visible
        %filterRange    % range of filter-data       
    end

   methods
       function  h = CellDataPatch(G, colorData, varargin)
           opt = struct('nbins',        100, ...
                        'colorData',     [], ...
                        'Parent',        [], ...
                        'fastMode',   false,  ...
                        'filterData',    [], ...
                        'playMode',  false, ...
                        'logScale',  false);  
                    
           if mod(numel(varargin), 2) == 1   % cell-subset is provided
               cells = varargin{1};
               if ~islogical(cells), cells = find(cells); end
               [opt, patchOpt] = merge_options(opt, varargin{2:end});
           else
               cells = (1:G.cells.num)';   % assume all cells
               [opt, patchOpt] = merge_options(opt, varargin{:});
           end
           
           if isempty(opt.Parent)
               opt.Parent = gca;
           end

           h.patchOpt = [{'Parent', opt.Parent}, patchOpt];
           
           if isempty(colorData)  % just use a vector of ones ...
               colorData = ones(G.cells.num,1);
           end
           
           if isempty(opt.filterData) % assume colordata and filterdata are the same
              opt.filterData = colorData;
           end
               
           assert(numel(colorData)==G.cells.num)
           assert(numel(opt.filterData)==G.cells.num)
            
           % store data in hidden struct to enable set/get corresponding properties
           h.gridData = struct('G', G, ...       % don't need all, but ...
                               'cells',      cells, ...
                               'colorData',  colorData, ...
                               'filterData', opt.filterData);
                                  
           h.fastMode = opt.fastMode;
           h.nbins    = opt.nbins;
           h.logScale = opt.logScale;
           % create empty patches with supplied options
           h = setEmptyMainPatch(h);
           %h.patchMain   = emptyPatch('Parent', opt.Parent, h.patchOpt{:});
           
           h = setEmptyAnimPatch(h);
%            if verLessThan('matlab', '8.4.0')
%                h.patchSlices = zeros(1, h.nbins);
%                h.patchLids   = zeros(1, h.nbins);
%            else
%                h.patchSlices = gobjects(1, h.nbins);
%                h.patchLids   = gobjects(1, h.nbins);
%            end
%            for k = 1:h.nbins
%                h.patchSlices(k) = emptyPatch('Parent', opt.Parent, h.patchOpt{:});
%                h.patchLids(k)   = emptyPatch('Parent', opt.Parent, h.patchOpt{:});
%            end
               
           if ~opt.playMode
               h.patchType = 'single';
           else
               h.patchType = 'multiple';
           end
           
           h.createPatch;
       end
           
       function createPatch(h)
           if ~h.playMode
               delete(h.patchSlices);
               delete(h.patchLids);
               h.main_ix = binPartitionPatchFaces(h.gridData.G, h.gridData.cells, ...
                                                  h.gridData.filterData, 1);  % single bin
                                              
               set(h.patchMain, 'Faces',          h.main_ix.faces, ...
                                'Vertices',       h.main_ix.vertices, ...
                                'FaceVertexCData', h.colorData(h.main_ix.faceCells,:), ...
                                'Visible', 'on');
           else % multiple
               assert(~isempty(h.gridData.filterData), 'FilterData empty, can''t make animation patch...\n')
               h.patchMain.Visible = 'off';
               h = setEmptyAnimPatch(h);
              [h.slice_ix, h.lid_ix, h.binMeans] = ...
                  binPartitionPatchFaces(h.gridData.G, h.gridData.cells, ...
                                         h.gridData.filterData, h.nbins, h.logScale, h.sliceThreshold);
              for k = 1:h.nbins
                  set(h.patchSlices(k), 'Faces',          h.slice_ix(k).faces, ...
                                        'Vertices',       h.slice_ix(k).vertices, ...
                                        'FaceVertexCData', h.colorData(h.slice_ix(k).faceCells,:), ...
                                        'Visible', 'off');   
                  set(h.patchLids(k), 'Faces',          h.lid_ix(k).faces, ...
                                      'Vertices',       h.lid_ix(k).vertices, ...
                                      'FaceVertexCData', h.colorData(h.lid_ix(k).faceCells,:), ...
                                      'Visible', 'off');
                                  % set only first to visible
              end
              h.patchSlices(1).Visible = 'on';
              h.patchLids(1).Visible = 'on';
           end
       end
       
       function step(h, n)  % step n bins (forward or backward)
           if h.playMode
               cb = h.curBin;
               n = max(min(n, h.nbins-cb), -cb);
               if n~=0
                   if n > 0
                       for k = (cb+1):(cb+n)
                           h.patchSlices(k).Visible = 'on';
                           h.patchLids(k).Visible = 'on';
                           if k-h.sliceThreshold > 0
                               h.patchLids(k-h.sliceThreshold).Visible = 'off';
                           end
                       end
                   else % n negative
                       for k = cb:-1:(cb+n+1)
                           h.patchSlices(k).Visible = 'off';
                           h.patchLids(k).Visible = 'off';
                           if k-h.sliceThreshold > 0
                               h.patchLids(k-h.sliceThreshold).Visible = 'on';
                           end
                       end
                   end
               end
           end
       end                   
               
       function set.cells(h, val)
           if islogical(val)
               val = find(val);
           end
           h.gridData.cells = val;
           h.createPatch;
       end
       function val = get.cells(h)
           val = h.gridData.cells;
       end
       
       function set.colorData(h, val)
           assert(size(val,1) == h.gridData.G.cells.num);
           h.gridData.colorData = val;
           if ~h.playMode  % single patch
               set(h.patchMain,'FaceVertexCData', h.colorData(h.main_ix.faceCells,:));
           else
               for k = 1:h.nbins
                   set(h.patchSlices(k), 'FaceVertexCData', h.colorData(h.slice_ix(k).faceCells,:)); 
                   set(h.patchLids(k),   'FaceVertexCData', h.colorData(h.lid_ix(k).faceCells,:));  
               end
           end
       end
       function val = get.colorData(h)
           val = h.gridData.colorData;
       end
       
       function set.filterData(h, val)
           assert(numel(val) == h.gridData.G.cells.num);
           h.gridData.filterData = val;
           h.createPatch;
       end
       function val = get.filterData(h)
           val = h.gridData.filterData;
       end
       
       function set.value(h, val)
           if ~h.playMode  % should have no effect
               warning('Setting the value has no effect in single patch mode')
           else
               cb     = h.curBin;
               bin_ix = find(h.binMeans <= val, 1, 'last');
               if isempty(bin_ix)
                   bin_ix = 0;
               end
               h.step(bin_ix-cb);
           end
       end
       function val = get.value(h)
           val = h.binMeans(h.curBin);
       end
       
       function set.curBin(h, val)
           if val > 0
               h.value = h.binMeans(val);
           else
               h.value = h.binMeans(1)-1;
           end
       end
       function val = get.curBin(h)
           str = {h.patchSlices(:).Visible};
           val = find(strcmp('on', str), 1, 'last');
           if isempty(val)
               val = 0;
           end
       end
       
       function  set.playMode(h, val)
           if val ~= h.playMode % playMode has changed
               if val  % switch  to playmode
                   h.patchType = 'multiple';
               else    % switch to single and clear multiple 
                   h.patchType = 'single';
                  % h = setEmptyAnimPatch(h);
               end
               h.createPatch;
           end
       end
      function val = get.playMode(h)
          val = strcmp(h.patchType, 'multiple');
      end
      
      function set.Visible(h, val)
          if h.playMode % set all
              if strcmp(val, 'on')
                  h.step(h.nbins);
              else
                  h.step(-h.nbins);
              end
          else
              h.patchMain.Visible = val;
          end
      end
      function val = get.Visible(h)
          if h.playMode
              if h.curBin == h.nbins
                  val = 'on';
              elseif h.curBin == 0
                  val = 'off';
              else
                  val = 'mixed';
              end
          else
              val = h.patchMain.Visible;
          end
      end
   end % methods
end % classdef

% -------------------------------------------------------------------------

% helpers
%function h = setEmpty(h)
%set(h, 'Vertices', [], 'Faces', []);
%end
function p = emptyPatch(opt)
p = patch('Vertices', [], 'Faces', [], 'FaceVertexCData', [], ...
          'FaceColor', 'flat', 'Visible', 'off', opt{:});
end

function h = setEmptyMainPatch(h)
h.patchMain = emptyPatch(h.patchOpt);
end

function h = setEmptyAnimPatch(h)
if verLessThan('matlab', '8.4.0')
    h.patchSlices = zeros(1, h.nbins);
    h.patchLids   = zeros(1, h.nbins);
else
    h.patchSlices = gobjects(1, h.nbins);
    h.patchLids   = gobjects(1, h.nbins);
end
for k = 1:h.nbins
    h.patchSlices(k) = emptyPatch(h.patchOpt);
    h.patchLids(k)   = emptyPatch(h.patchOpt);
end
end

function [slice_ix, lid_ix, binCenters] = binPartitionPatchFaces(G, cells, data, nbins, logScale, sliceThreshold)
if nargin < 6
    sliceThreshold = round(sqrt(nbins)); % balance between storage and max speed
end
assert(numel(data) == G.cells.num);
slice_ix = repmat(struct('faces', [], 'vertices', [], 'faceCells', []), nbins, 1);
lid_ix   = repmat(struct('faces', [], 'vertices', [], 'faceCells', []), nbins, 1);

[nodes, pos] = deal(G.faces.nodes, G.faces.nodePos);
if nbins == 1 % single patch, just one lid
    if numel(cells) > 0
        [f, fc] = boundaryFaces(G, cells);
        [slice_ix.faces, vertix] = get_face_topo(nodes, pos, f);
        slice_ix.vertices = G.nodes.coords(vertix, :);
        slice_ix.faceCells = fc;
        binCenters = []; % don't need this
    end
else % more than one bin
    maxn = double(max(nodes));
    data_local = data(cells);
    
    [minVal, maxVal] = deal(min(data_local), max(data_local));
    % make some slack in case degenerate case
    dv = max( max(abs([maxVal, minVal]))*sqrt(eps), sqrt(eps));
    [minVal, maxVal] = deal(minVal-dv, maxVal+dv);
    
    if ~logScale
        bins = minVal + (0:nbins)*(maxVal-minVal)/nbins;
        binCenters = (bins(1:end-1)+bins(2:end))/2;
    else
        logbins = log10(minVal) +  (0:nbins)*(log10(maxVal)-log10(minVal))/nbins;
        bins    = 10.^logbins;
        binCenters = 10.^((logbins(1:end-1)+logbins(2:end))/2);
    end
    
    p = zeros(G.cells.num, 1); % partitioin vector
    %p(cells) = ceil(nbins*(data_local-minVal)/(maxVal-minVal));
    [~, ~, p(cells)] = histcounts(data_local, bins);
    CG = generateCoarseGrid(G, p+1);
    p0 = [0; p];
    [nc, connPos] = deal(CG.faces.neighbors, CG.faces.connPos);
    connPos = [connPos(1:end-1), connPos(2:end)-1];
    % treat exterior and grid outside region as the same, set to 1
    nc = max(nc,1);
    nc = nc-1;  % reset to correct number
    swap =  or(nc(:,1)>nc(:,2), nc(:,1) == 0);
    nc(swap, :) = nc(swap, [2 1]);
    % sort
    [nc, ix] = sortrows(nc);
    connPos  = connPos(ix,:);

    % gather
    [cc, nf] = rlencode(nc(:,1));
    if cc(1) == 0   % remove extra faces from outside cells
        nc = nc(nf(1)+1:end,:);
        connPos = connPos(nf(1)+1:end,:);
        [cc, nf] = deal(cc(2:end), nf(2:end));
    end
    if numel(cc) < nbins   % include empty bins
        nf_tmp = nf;
        nf     = zeros(nbins,1);
        nf(cc) = nf_tmp;
    end
    fPos = [1; cumsum(nf)+1];
    % now loop through
    for k  = 1:nbins
        cix  = fPos(k):(fPos(k+1)-1);
        curn = nc(cix, 2);
        curPos = connPos(cix,:);

        % select indices to fine faces
%         f = cellfun(@(x)CG.faces.fconn(x(1):x(2)), ...
%                     mat2cell(curPos, ones(size(curPos,1),1), 2), ...
%                     'UniformOutput', false);
%                 
%         % ...
        mm =  mat2cell(curPos, ones(size(curPos,1),1), 2);
        f = cellfun(@(x)CG.faces.fconn(x(1):x(2)), mm, ...
                    'UniformOutput', false);        
                
        % pick which become slices
        %nslice = sum((curn-k) < sliceThreshold);
        slix = or((curn-k) > sliceThreshold, curn == 0);
%        fs = vertcat(f{1:nslice});
        fs = vertcat(f{slix});
        [slice_ix(k).faces, vertix] = get_face_topo(nodes, pos, fs, maxn);
        slice_ix(k).vertices = G.nodes.coords(vertix, :);
        slice_ix(k).faceCells = getFaceCells(G, fs, p0, k);
        
        %fl = vertcat(f{nslice+1:end});
        fl = vertcat(f{~slix});
        [lid_ix(k).faces, vertix] = get_face_topo(nodes, pos, fl, maxn);
        lid_ix(k).vertices = G.nodes.coords(vertix, :);
        lid_ix(k).faceCells = getFaceCells(G, fl, p0, k);
    end
end
end

function c = getFaceCells(G, f, p0, cc)
n  = G.faces.neighbors(f,:);
c  = n(:,1);
ii = p0(n(:,2)+1) == cc;
c(ii) = n(ii,2);
end

function [f, present] = get_face_topo(nodes, pos, faces, maxn)
if ~isempty(faces)
    if nargin < 4
        maxn =  double(max(nodes));
    end
   eIX = pos;
   nn  = double(diff([pos(faces), ...
                      pos(faces + 1)], [], 2));
   fn  = double(nodes(mcolon(eIX(faces), eIX(faces + 1) - 1), 1));

   m   = numel(faces);
   n   = max(nn);
   f   = nan([n, m]);

   % Extract only those nodes/vertices actually present in the subset of
   % grid faces represented by 'faces'.  Create local numbering for these
   % vertices.
   %

   %present           = false([maxn, 1]);
   %present(fn)       = true;

   present = sparse(fn, 1, 1, maxn, 1) > 0;
   %pp = pp > 0;
   %node_num          = zeros([maxn, 1]);
   %node_num(present) = 1 : sum(double(present));


   node_num = sparse(find(present), 1, (1:nnz(present))',  maxn, 1);
   off = reshape((0 : m - 1) .* n, [], 1);

   f(mcolon(off + 1, off + nn)) = node_num(fn);

   %tmp         = isfinite(f);
   %nnode       = sum(tmp,1);
   %ind         = sub2ind(size(f),nnode,1:size(f,2));
   %tmp         = repmat(f(ind),size(f,1),1);
   %f(isnan(f)) = tmp(isnan(f));
   % PATCH requires that the 'Faces' property be a matrix of size
   % (number of faces)-by-(number of vertices).
   %
   f = f .';
else
    f = [];
    present = [];
end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
