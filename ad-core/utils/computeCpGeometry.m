function [cellCenters, cellFaceCenters, cellDims]=computeCpGeometry(G,grdecl,varargin)
opt=struct('on_all_cpc',false);
opt=merge_options(opt,varargin{:});
[x, y, z] = buildCornerPtNodes(grdecl);

px1 = x(1:2:end, 1:2:end, 1:2:end);
px2 = x(2:2:end, 1:2:end, 1:2:end);
px3 = x(1:2:end, 2:2:end, 1:2:end);
px4 = x(2:2:end, 2:2:end, 1:2:end);
px5 = x(1:2:end, 1:2:end, 2:2:end);
px6 = x(2:2:end, 1:2:end, 2:2:end);
px7 = x(1:2:end, 2:2:end, 2:2:end);
px8 = x(2:2:end, 2:2:end, 2:2:end);
clear x

py1 = y(1:2:end, 1:2:end, 1:2:end);
py2 = y(2:2:end, 1:2:end, 1:2:end);
py3 = y(1:2:end, 2:2:end, 1:2:end);
py4 = y(2:2:end, 2:2:end, 1:2:end);
py5 = y(1:2:end, 1:2:end, 2:2:end);
py6 = y(2:2:end, 1:2:end, 2:2:end);
py7 = y(1:2:end, 2:2:end, 2:2:end);
py8 = y(2:2:end, 2:2:end, 2:2:end);
clear y

pz1 = z(1:2:end, 1:2:end, 1:2:end);
pz2 = z(2:2:end, 1:2:end, 1:2:end);
pz3 = z(1:2:end, 2:2:end, 1:2:end);
pz4 = z(2:2:end, 2:2:end, 1:2:end);
pz5 = z(1:2:end, 1:2:end, 2:2:end);
pz6 = z(2:2:end, 1:2:end, 2:2:end);
pz7 = z(1:2:end, 2:2:end, 2:2:end);
pz8 = z(2:2:end, 2:2:end, 2:2:end);
clear z

p1 = [px1(:) py1(:) pz1(:)];
clear px1 py1 pz1
p2 = [px2(:) py2(:) pz2(:)];
clear px2 py2 pz2
p3 = [px3(:) py3(:) pz3(:)];
clear px3 py3 pz3
p4 = [px4(:) py4(:) pz4(:)];
clear px4 py4 pz4
p5 = [px5(:) py5(:) pz5(:)];
clear px5 py5 pz5
p6 = [px6(:) py6(:) pz6(:)];
clear px6 py6 pz6
p7 = [px7(:) py7(:) pz7(:)];
clear px7 py7 pz7
p8 = [px8(:) py8(:) pz8(:)];
clear px8 py8 pz8

xfc1 = (p1 + p3 + p5 + p7)/4;
xfc2 = (p2 + p4 + p6 + p8)/4;
yfc1 = (p1 + p2 + p5 + p6)/4;
yfc2 = (p3 + p4 + p7 + p8)/4;
zfc1 = (p1 + p2 + p3 + p4)/4;
zfc2 = (p5 + p6 + p7 + p8)/4;
clear p1 p2 p3 p4 p5 p6 p7 p8

cc = (zfc1 + zfc2)/2;

if ~(opt.on_all_cpc)
   cellNo     = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   cartCellNo = G.cells.indexMap(cellNo);
   fType      = G.cells.faces(:, 2);
   cellCenters      =  cc( G.cells.indexMap, :);
else
   nc=prod(grdecl.cartDims);
   cartCellNo = rldecode(1:nc, 6, 2) .';
   fType =repmat([1:6]',nc,1);
   cellCenters=cc;
end


cellFaceCenters  =  ...
   ((fType == 1)*[1 1 1]).*xfc1(cartCellNo, :) + ...
   ((fType == 2)*[1 1 1]).*xfc2(cartCellNo, :) + ...
   ((fType == 3)*[1 1 1]).*yfc1(cartCellNo, :) + ...
   ((fType == 4)*[1 1 1]).*yfc2(cartCellNo, :) + ...
   ((fType == 5)*[1 1 1]).*zfc1(cartCellNo, :) + ...
   ((fType == 6)*[1 1 1]).*zfc2(cartCellNo, :);

if nargout > 2
    vnorm = create_vector_norm_function();
    cellDims = [vnorm(xfc2-xfc1), vnorm(yfc2-yfc1), vnorm(zfc2-zfc1)];
end
end

%--------------------------------------------------------------------------

function vnorm = create_vector_norm_function()
    if ~exist('verLessThan', 'file') || verLessThan('matlab', '9.3.0')
        % Fall-back for 'vecnorm'; introduced in MATLAB 9.3.0 (R2017b)
        vnorm = @(x) sqrt(sum(x.^2, 2));
    else
        vnorm = @(x) vecnorm(x, 2, 2);
    end
end
