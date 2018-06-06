function fgnum=findFracGridnum(G,cellnum)
% Find FracGrid number given global cell number

fgnum=1;
numfgrid=length(fieldnames(G.FracGrid));

fcellstart=G.FracGrid.Frac1.cells.start;
while cellnum>=fcellstart
    fgnum=fgnum+1;
    if fgnum>numfgrid
        break;
    end
    fcellstart=G.FracGrid.(['Frac',num2str(fgnum)]).cells.start;
end
fgnum=fgnum-1;

end