function reg = handleRegions(deck, G, varargin)

if ~isempty(G)
    an = G.cells.indexMap;    
elseif isfield(deck.GRID, 'ACTNUM')
    an = find(deck.GRID.ACTNUM);
else
    an = ':';
end
  
  
% also need actnum/G.cells.indexmap
% PVT-REGIONS
ntpvt = 1;
if isfield(deck.REGIONS, 'PVTNUM')
    reg.PVTNUM = deck.REGIONS.PVTNUM(an);
    ntpvt = max(reg.PVTNUM);
end
if ntpvt == 1
    reg.PVTNUM = [];
    reg.PVTINX = ':';
else
    reg.PVTINX = cellfun(@(x)(x==reg.PVTNUM), num2cell(1:ntpvt), 'UniformOutput', false);
end

% SAT-REGIONS AND POSSIBLY SURF-REGIONS
one_region = true;
if isfield(deck.REGIONS, 'SATNUM')
    reg.SATNUM = deck.REGIONS.SATNUM(an);
    if isfield(deck.REGIONS, 'SURFNUM')
       reg.SURFNUM = deck.REGIONS.SURFNUM(an);
       ntsatsurfact = max(max(reg.SATNUM), max(reg.SURFNUM));
       reg.SATINX = cellfun(@(x)find(x==reg.SATNUM), num2cell(1:ntsatsurfact), 'UniformOutput', ...
                            false);
       reg.SURFINX = cellfun(@(x)find(x==reg.SURFNUM), num2cell(1:ntsatsurfact), 'UniformOutput', ...
                             false);
       one_region = false;
    else
       ntsat = max(reg.SATNUM);
       if ntsat > 1
          reg.SATINX = cellfun(@(x)find(x==reg.SATNUM), num2cell(1:ntsat), 'UniformOutput', ...
                               false);
          one_region = false;
       end
    end
elseif isfield(deck.REGIONS, 'SURFNUM')
      error('SATNUM keyword required when surfactant is used.');
end

if one_region
      reg.SATNUM = [];
      reg.SATINX = ':';
end

% IMB-REGIONS
ntsat = 1;
if isfield(deck.REGIONS, 'IMBNUM')
    reg.IMBNUM = deck.REGIONS.IMBNUM(an);
    ntsat = max(reg.IMBNUM)/2;
end
if ntsat == 1
    reg.IMBNUM = [];
    reg.IMBINX = ':';
else
    reg.IMBINX = cellfun(@(x)find(x==reg.IMBNUM), num2cell(1:ntsat), 'UniformOutput', false);
end

% ROCK-REGIONS
ntrocc = 1;
if isfield(deck.REGIONS, 'ROCKNUM')
    reg.ROCKNUM = deck.REGIONS.ROCKNUM(an);
    ntrocc = max(reg.ROCKNUM);
end
if ntrocc == 1
    reg.ROCKNUM = [];
    reg.ROCKINX = ':';
else
    reg.ROCKINX = cellfun(@(x)find(x==reg.ROCKNUM), num2cell(1:ntrocc), 'UniformOutput', false);
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
