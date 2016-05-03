function ht = my_yticklabels(varargin)

% ht = my_yticklabels(Ha, xtickpos, xtickstring)
% or
% ht = my_yticklabels(xtickpos, xtickstring)

% Pekka Kumpulainen 12.2.2008

textopts = {};
if length(varargin{1})==1 && ...
        ishandle(varargin{1}) && ...
        strcmpi(get(varargin{1},'Type'),'axes');
    Ha = varargin{1};
    ytickpos = varargin{2};
    ytickstring = varargin{3};
    if nargin > 3
        textopts = varargin(4:end);
    end
else
    Ha = gca;
    Hfig = get(Ha,'Parent');
    ytickpos = varargin{1};
    ytickstring = varargin{2};
    addy = varargin{3};
    addx = varargin{4};
    if nargin > 4
        textopts = varargin(3:end);
    end
end

%% Make YTickLabels
NTick = length(ytickpos);
Xbot = min(get(gca,'XLim'))+addx;
ytickpos = ytickpos + addy;
ht = zeros(NTick,1);
for ii = 1:NTick
    ht(ii) = text('String',ytickstring{ii}, ...
        'Units','data', ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'center ', ...
        'Position',[Xbot ytickpos(ii)], ...
        'Tag','MUXTL');
end
if ~isempty(textopts)
    set(ht,textopts{:})
end

%% squeeze axis if needed

set(Ha,'Units','pixels')
Axpos = get(Ha,'Position');
% set(Hfig,'Units','pixels')
% Figpos = get(Hfig,'Position');

set(ht,'Units','pixels')
TickExt = zeros(NTick,4);
for ii = 1:NTick
    TickExt(ii,:) = get(ht(ii),'Extent');
end

needmove = -(Axpos(1) + min(TickExt(:,1)));

if needmove>0;
    Axpos(1) = Axpos(1)+needmove+2;
    Axpos(3) = Axpos(3)-needmove+2;
    set(Ha,'Position',Axpos);
end

set(Ha,'Units','normalized')
set(ht,'Units','normalized') 