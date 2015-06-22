function checkIfCoplanar(fracplanes)

warning('off','MATLAB:triangulation:EmptyTri3DWarnId');
for i = 1:numel(fracplanes)
    x = [];
    y = [];
    z = [];
    for j = 1:size(fracplanes(i).points,1)
        x = [x;fracplanes(i).points(j,1)]; %#ok
        y = [y;fracplanes(i).points(j,2)]; %#ok
        z = [z;fracplanes(i).points(j,3)]; %#ok
    end
    assert(all(iscoplanar([x,y,z])),['Input Points must be coplanar. Point set ',...
        num2str(i),' is not!']);
end

warning('on','MATLAB:triangulation:EmptyTri3DWarnId');
return