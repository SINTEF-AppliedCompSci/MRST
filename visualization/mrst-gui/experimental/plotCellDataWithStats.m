function varargout = plotCellDataWithStats(G, data, varargin)

cells = (1 : G.cells.num) .';
if mod(numel(varargin), 2) == 1,
   if isnumeric(varargin{1}),
      cells = varargin{1};
   elseif islogical(varargin{1}) && ...
         numel(varargin{1}) == G.cells.num,
      cells = find(varargin{1});
   else
      error(['Third parameter ''cells'' must either be a list of ', ...
             'explicit cell indices or a logical mask into the '  , ...
             'grid''s cells.']);
   end
   varargin = varargin(2 : end);
end

h = plotCellData(G, data, cells, varargin{:});
pltax = gca;
set(pltax,'Position',[.05 .05 .55 .9]);

axes('Position',[.675 .1 .3 .4]);
[nv,bv] = hist(data,50);
hb = bar(bv,nv,1);
set(get(hb,'Children'),'CData',bv,'EdgeColor','none'); axis tight

S = sort(data);
N = size(data,1);
Q = interp1q([0 (0.5:(N-0.5))./N 1]',S([1 1:N N],:),[.25; .75]); 

annotation('textbox',[.675 .525 .3 .25], 'String', {...
    sprintf('%s: %d cells', inputname(2), numel(cells)), ...
    sprintf('Range  : [%.3g, %.3g]',min(data), max(data)), ...
    sprintf('Mean   : %g', mean(data)), ...
    sprintf('Std.dev: %g', std(data)), ...
    sprintf('Median : %g', median(data)), ...
    sprintf('Quartiles: %.3g, %.3g', Q), ...
    });
set(gcf, 'CurrentAxes', pltax);
if nargout > 0, varargout{1} = h; end