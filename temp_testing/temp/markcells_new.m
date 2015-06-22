function [fraCells, remove] = markcells_new(G,GlobTri,fracplanes)
remove = [];
Gtri = GlobTri.Tri;
map = GlobTri.map;
fraCells = cell(numel(fracplanes),1);
possible_ratios = 0:0.0001:1;
possible_ratios = possible_ratios(mod(1,possible_ratios)==0);
count = 0;
for i = 1:numel(fracplanes)
    count = count+1;
    ratio = max(G.faces.areas)/(10*polyArea3D(fracplanes(i).points));
    [~,loc] = min(abs(possible_ratios-ratio));
    ratio = possible_ratios(loc);
    if ratio == 1
        warning(['Fracture ',num2str(i),' is not a long fracture and ',...
            'will be removed from further calculations.']);
        remove = [remove;i]; %#ok
        count = count-1;
        continue
    end
    [p,in,on] = pointsInPolygonGenerator(fracplanes(i).points,'ratio',ratio);
    cells = unique(map(pointLocation(Gtri,p(in~=on,:))));
    fraCells{count,1} = cells;
end
return



