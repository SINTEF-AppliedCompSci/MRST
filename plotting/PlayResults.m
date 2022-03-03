function PlayResults(fig)
    figTitle = fig.title;
    figStyle = fig.style;
    figTag   = fig.tag;
    states   = fig.states;
    G        = fig.grid;  
    h = findobj( 'Type', 'Figure', 'Name', figTitle);
    if(isempty(h))
        figure('Name',        figTitle, ...
               'Tag',         figTag,   ...
               'NumberTitle', 'off',     ...
               'WindowStyle', figStyle );
    else
        figure(h);   
    end                 
    plotToolbar(G, states, 'field', 's:1', 'plot1d', true, ...
                'lockCaxis', true, 'startplayback', true);
end