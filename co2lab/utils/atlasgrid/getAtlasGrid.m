function varargout = getAtlasGrid(varargin)
% Get GRDECL grids and datasets for CO2 Atlas datasets
%
% SYNOPSIS:
%   getAtlasGrid();    % lists all possible datasets
%
%   [grdecls datasets petroinfo] = getAtlasGrid();
%   [grdecls datasets petroinfo] = getAtlasGrid({'Johansenfm', 'Brynefm'});
%   [grdecls datasets petroinfo] = getAtlasGrid('Johansenfm');
%
%   [grdecls ...] = getAtlasGrid(pn1, pv1, ..)
%   [grdecls ...] = getAtlasGrid('Johansenfm', pn1, pv1, ...)
%   [grdecls ...] = getAtlasGrid({'Johansenfm', ..}, pn1, pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of optional property names/property values:
%                   
%    - coarsening: Coarsening factor. If set to one, a grid with approximately
%                  one cell per datapoint is produced. If set to two, every second
%                  datapoint in x and y direction is used, giving a reduction to
%                  1/4th size. Default: 1
%     - refining:  Refining factor. If set to an integer n > 1, each 
%                  original cell is laterally split up into n x n cells. 
%
%    - nz: Vertical / k-direction number of layers. Default: 1
%
%    - make_deck: Boolean. If set to true, the routine will create cell
%             arrays of data members in ECLIPSE grid structures to
%             represent the 3D grid. If false, the routine will create a
%             structure that contains node points, depth of top surface,
%             thickness of formation, etc. 
%
% NOTE:
%     If called without parameters, the function will list all files in the
%     data directory.
%
% RETURNS:
%   grdecls - Cell arrays of data members in ECLIPSE grid structures
%             suitable for processGRDECL if 'make_deck' is true.
%             Otherwise, a data structure that contains arrays giving the
%             nodes, depth of top surface, thickness of formation, etc.
%
%   datasets - (OPTIONAL) raw datasets used to produce GRDECL structs.
%
%   petroinfo - (OPTIONAL) average perm / porosity as given in the atlas.
%             Will contain NaN where values are not provided or possible
%             to calculate.
%
% SEE ALSO:
%   `convertAtlasTo3D`, `downloadDataSets`

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
if nargout == 0
   n = getNames();
   n = n( cellfun(@(x) ~isempty(strfind(x,'_')), n) );
   n = regexp(n, '_', 'split');
   n = unique(cellfun(@(x) x{1}, n, 'UniformOutput', false));
   disp 'Available grids are: '
   fprintf('\t- %s\n', n{:});
   return
end

opt = struct('nz',              1           , ...
             'coarsening',      1           , ...
             'refining',        1           , ...
             'Verbose',         mrstVerbose , ...
             'make_deck',       true);
if mod(nargin, 2) == 1
   names = varargin{1};
   if nargin > 1
      varargin = {varargin{2:end}}; %#ok
   else
      varargin = {};
   end
else
   names = {};
end

opt = merge_options(opt, varargin{:});

if ischar(names)
   names = {names};
end


% Handle grids belonging to the Jurassic formations
% This is done by sampling the Jurassic formation for the top heights
jurassic_mid = {'Brentgrp',...
   'Brynefm',...
   'Sleipnerfm'};

jurassic_top = {'Sognefjordfm',...
   'Fensfjordfm',...
   'Krossfjordfm',...
   'Huginfmeast',...
   'Huginfmwest',...
   'Sandnesfm',...
   'Ulafm'};

if ~isempty(names)
   juractive = [0 0];
   
   jurassic_mid = intersect(jurassic_mid, names);
   jurassic_top = intersect(jurassic_top, names);
   if ~isempty(jurassic_mid)
      names = [names, {'Jurassicmid'}];
      juractive(1) = 1;
   end
   if ~isempty(jurassic_top)
      names = [names, {'Jurassic'}];
      juractive(2) = 1;
   end
else
   juractive = [1 1];
end

datasets = readAtlasGrids(names, opt.coarsening);

jura_top = {jurassic_mid, jurassic_top};
jurassic = {'Jurassicmid', 'Jurassic'};
if(opt.refining>1)
   for i=1:numel(datasets)
      datasets{i}=refineInfo(datasets{i},opt.refining);
   end
end
findname = @(name) datasets{cellfun(@(x) strcmpi(x.name, name), datasets)};
   
[decks, petroinfo] = deal({});

for i = 1:2
   jur = jura_top{i};
   if ~juractive(i)
      continue
   end
   jtop = findname(jurassic{i});
   
   for j = 1:numel(jur)
      jthick = findname(jur{j});
      dispif(opt.Verbose, ['Processing ' jur{j} ' ...\n']);
      if(opt.make_deck)
         deck = convertAtlasTo3D(jthick.meta, jtop.meta, jthick.data, jtop.data, opt.nz);
      else
         deck = convertAtlasToStruct(jthick.meta, jtop.meta, jthick.data, jtop.data);
      end
      deck.name = jthick.name;
      
      petroinfo{end+1} = getAvgRock(deck.name); %#ok
      % NB: if heterogeneous rock property data becomes available for any
      % formation belonging to jurassic, call to updateWithHeterogenVals()
      % should be inserted here.
      decks{end+1} = deck;%#ok
   end
end

% Handle grids with seperate thickness / top maps of same name
for i = 1:numel(datasets)
   top = datasets{i};
   if ~strcmpi(top.variant, 'top')
      continue
   end
   
   for j = 1:numel(datasets)
      thick = datasets{j};
      if strcmpi(thick.variant, 'thickness') && strcmpi(thick.name, top.name)
         dispif(opt.Verbose, ['Processing ', top.name, ' ...\n']);
         if(opt.make_deck)
            deck = convertAtlasTo3D(thick.meta, top.meta, thick.data, top.data, opt.nz);
         else
            deck = convertAtlasToStruct(thick.meta, top.meta, thick.data, top.data);
         end
         deck.name = thick.name;
         
         % first get average rock properties as reported in NPD's Atlas
         petroinfo{end+1} = getAvgRock(deck.name); %#ok
         
         % then get heterogeneous rock properties (if they exist).
         % They are added to deck structure, and average rock properties in
         % petroinfo structure are updated.
         [deck, petroinfo{end}] = updateWithHeterogeneity(deck, ...
                                  petroinfo{end}, top.meta, opt);
         
         decks{end+1} = deck; %#ok
         break;
      end
   end
end
varargout{1} = decks;
if nargout > 1
   varargout{2} = datasets;
   if nargout > 2
      varargout{3} = petroinfo;
   end
end
end

function target_path = select_filepath(gdirs, gname)
   assert (ischar(gname), 'Internal Error');

   for p = reshape(gdirs, 1, [])

      files = dir(p{1});
      if any(strcmpi({files.name}, gname))
         target_path = fullfile(p{1}, gname);
         return;
      end
   end
   error(['Could not find ' gname ' in any directory.']);
end

function datasets = readAtlasGrids(names, coarsening)

[grids, gdir] = getNames();

datasets = {};
for i = 1:numel(grids)
    g = grids{i};
    
    nstr = regexp(g, '_', 'split');

    % Designate the type based on the file names
    ind = strcmpi(nstr, 'thickness') | strcmpi(nstr, 'top');
    tmp.name = [nstr{~ind}];
    
    if any(strcmpi(nstr, 'thickness'))
        tmp.variant = 'thickness';
        isTop = false;
    else
        tmp.variant = 'top';
        isTop = true;
    end
    
    
    if ~isempty(names) && ~any(strcmpi(tmp.name, names))
        continue;
    end
    
    % Skip projections and mat files
    if any(strcmpi(g(end-3:end), {'.mat', '.prj'}))
        continue
    end
    % Skip any readme files
    if any(strcmpi(g, {'readme.txt', 'readme'}))
        continue
    end
    
    %[meta, data] = readAAIGrid(fullfile(gdir, g));
    [meta, data] = readAAIGrid(select_filepath(gdir, g));
    
    % The rawdata files of North Sea formations have neg depth values for
    % top surface, while the rawdata files of Barents Sea formations have
    % pos depth values for top surface. We want data to contain pos depth
    % values. Here, sign of data is checked and is switched if neg:
    if isTop && any(find(data<0))
        data = -data;
    end
    
    meta.cellsize = meta.cellsize*coarsening;
    
    tmp.meta = meta;
    tmp.data = data(1:coarsening:end, 1:coarsening:end);
    
    tmp.meta.ncols = size(tmp.data, 2);
    tmp.meta.nrows = size(tmp.data, 1);
    
    tmp.meta.dims = [tmp.meta.nrows, tmp.meta.ncols];
    
    datasets{end+1} = tmp; %#ok
    
end
end

function [n, gdir] = getNames()
   [n1, gdir1] = getNamesFromDataset('co2atlas');
   [n2, gdir2] = getNamesFromDataset('co2atlasnorwegiansea');
   [n3, gdir3] = getNamesFromDataset('co2atlasbarentssea');
   
   n = [n1, n2, n3];
   gdir = {gdir1, gdir2, gdir3};
end


function [n, gdir] = getNamesFromDataset(dataset)
    gdir = getDatasetPath(dataset);

    dir_grid = dir(gdir);
    n = {dir_grid(~[dir_grid.isdir]).name};
end


function avgrock = getAvgRock(name)
    if strcmpi(name(end-1:end), 'fm')
        name = name(1:end-2);
    end
    switch lower(name)
        
        % North Sea
        case 'utsira'
            tmp = [1000, 0.2112];
        case 'skade'
            tmp = [1000, 0.2112];
        case 'sandnes'
            tmp = [150, 0.0875];
        case {'sognefjord', 'fensfjord', 'krossfjord'}
            % chp 4, p. 62 of Compiled Atlas says: Sognefjord Delta
            % comprised of these three formations, partly separated by thin
            % shales, and treated as one aquifer unit
            tmp = [300 0.1949];
        case 'statfjord'
            tmp = [200 0.1071];
        case 'gassum'
            tmp = [450	0.1165];
        case 'farsund'
            tmp = [150	0.0960];
        case {'johansen', 'cook'}
            % p. 29: porosities and permeabilities in the order of 15-24%
            % and 100-1000 mD respectively have been recorded
            tmp = [300 0.20];
        case 'fiskebank'
            tmp = [1000	0.2500];
        case 'stord'
            tmp = [15 0.0600];
        case {'huginfmeast', 'huginfmwest'}
            % chp 4, p. 44 of Compiled Atlas says: average porosity and
            % permeability is between 16-20% and 0.1-4000 mD respectively.
            % Only values for Hugin East are provided (see chp 4, p. 72).
            tmp = [500 0.1254]; % taken from North Sea Atlas, p. 62
        case 'bryne'
            tmp = [150	0.1280];
        case 'brentgrp'
            % chp 4, p. 38 of Compiled Atlas: this group is located at a
            % wide range of depths, so there is a complex distribution of
            % porosity and permeability. Exact values not stated in Atlas.
            tmp = [NaN NaN];
        case 'ula'
            % chp 4, p. 47 of Compiled Atlas says: porosities and
            % permeabilities are reported in the range of 15-22% and
            % 0.2-2800 mD respectively
            tmp = [300 0.185];
        case 'pliocenesand'
            % Exact values not stated in Atlas, but this formation is said
            % to be in communication with Skade (chp 4, p. 52)
            tmp = [NaN NaN];
            
        % Barents Sea (chp 6 of "CO2 Storage Atlas: Norwegian
        % Contential Shelf" from NPD)
        case 'bjarmeland'
            tmp = [300 0.23];   % p. 130: permeabilities in the order of
                                % 5-1000 mD have been recorded
                                % NB: Heterogeneous perm, poro, ntg data
                                % available for Prospect A
        case 'sto'
            tmp = [500 0.15];   % pg 128
        case 'nordmela'
            tmp = [1 0.15];     % pg 128
        case 'tubaen'
            tmp = [500 0.15];   % pg 128
            
        % Norwegian Sea (chp 5 of "CO2 Storage Atlas: Norwegian
        % Contential Shelf" from NPD)
        case 'tilje'
            tmp = [140 0.21 0.30]; % pg. 102 (including net-to-gross)
        case 'are'
            tmp = [140 0.21 0.30]; % pg. 102 (including net-to-gross)
        case 'garn'
            tmp = [580 0.27 0.25]; % pg. 103 (including net-to-gross)
        case 'ile'
            tmp = [580 0.27 0.25]; % pg. 103 (including net-to-gross)
        case 'not'
            tmp = [NaN NaN]; % pg. 87, 96: Not is a sealing formation
        case 'ror'
            tmp = [NaN NaN]; % pg. 84, 96: Ror is a sealing formation

        otherwise
            tmp = [NaN NaN];
    end
    avgrock.avgperm = tmp(1)*milli*darcy;
    avgrock.avgporo = tmp(2);
    if numel(tmp) == 3
       avgrock.avgntg = tmp(3);
    end
end
