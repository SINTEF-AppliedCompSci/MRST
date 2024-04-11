function set_fig_standardformat(fig, fig_title)

if ~iscell(fig_title)
    fig_title = {fig_title};
end

plots = get(fig, 'children');
set(gcf, 'color', 'white');
N = numel(plots);
for i = 1:N
    axis tight;
    pix = N-i+1;
    set(plots(pix), 'fontsize', 14);
    
    set(get(plots(pix), 'title'), 'String',  fig_title{i});
    set(get(plots(pix), 'xlabel'), 'String', 'meter');
    set(get(plots(pix), 'ylabel'), 'String', 'meter');
    set(get(plots(pix), 'zlabel'), 'String', 'meter');
end


if N == 1
    set(fig, 'position', [420, 70, 1200, 360]);
elseif N == 2
    set(fig, 'position', [300, 300, 1200, 400*numel(plots)])
else
    % multiple plots, we need larger figure
    set(fig, 'position', [300, 300, 2000, 800])
end

