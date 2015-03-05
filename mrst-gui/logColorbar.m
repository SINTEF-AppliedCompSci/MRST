function h = logColorbar(varargin)
    h = colorbar(varargin{:});
    
    ticks = get(h, 'YTick');
    
    newticks = arrayfun(@(x) ['1e', num2str(x)], ticks, 'UniformOutput', false);
    
    set(h, 'YTickLabel', newticks)
end
