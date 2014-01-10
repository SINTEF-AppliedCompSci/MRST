function reg = handleRegions(deck, G)
if nargin == 2
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
    reg.PVTINX = cellfun(@(x)find(x==reg.PVTNUM), num2cell(1:ntpvt), 'UniformOutput', false);
end

% SAT-REGIONS
ntsat = 1;
if isfield(deck.REGIONS, 'SATNUM')
    reg.SATNUM = deck.REGIONS.SATNUM(an);
    ntsat = max(reg.SATNUM);
end
if ntsat == 1
    reg.SATNUM = [];
    reg.SATINX = ':';
else
    reg.SATINX = cellfun(@(x)find(x==reg.SATNUM), num2cell(1:ntsat), 'UniformOutput', false);
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



