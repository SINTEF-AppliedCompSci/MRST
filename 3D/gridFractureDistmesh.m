function Gf = gridFractureDistmesh(G, fracplane, scaledplane, h, varargin)

opt = struct('type'         ,  'triangle', ...
             'useDistmesh'  ,   true     , ...
             'rectangular'  ,   false    , ...
             'minTriangles' ,   10       , ...
             'scale'        ,   2        );
opt = merge_options(opt, varargin{:});
%
flag = 1*strcmp(opt.type,'triangle') + 2*strcmp(opt.type,'pebi');
%
xyp = rotatePlane(scaledplane.points,[0 0 1]);
z = xyp(1,3);
xyp = xyp(:,1:2);
h = (h<1)*h^2 + (h>1)*sqrt(h); % assuming esize<0
h = h/opt.scale;
hstr = num2str(h,'%f');
sigfigs = hstr(hstr~='0' & hstr~='.');
nsigfigs = length(sigfigs); if nsigfigs>5, nsigfigs = nsigfigs-2; end
h = round(h,nsigfigs);
if h<0.05, h = 0.05; end
expand = h*10;
bbox = [min(xyp(:,1)) - expand, min(xyp(:,2)) - expand; ...
    max(xyp(:,1)) + expand, max(xyp(:,2)) + expand];

%-----------------------------triangle------------------------------------%
count = 0;
iter = 100;
while true
    count = count + 1;
    assert(count<10,' ');
    try
        if opt.rectangular
            fd = @ (p) drectangle(p, min(xyp(:,1)), max(xyp(:,1)), min(xyp(:,2)), max(xyp(:,2)));
            [p,t] = distmesh_2d(fd, @huniform, h, bbox, iter, ...
                xyp);
        else
            [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, ...
                xyp,xyp);
        end
        assert(size(t,1)>opt.minTriangles,'');
        close gcf
        Gf = triangleGrid(p,t);
        computeGeometry(Gf);
        computeGeometry(makeLayeredGrid(Gf,1));
        break;
    catch
        h = h/2;
        iter = iter+50;  
    end
end
%---------------------------------pebi------------------------------------%
iter = 50;
if flag == 2
    %     [Px, Py] = lloydsAlgorithm(p(:,1),p(:,2), xyp, 200, true);
    while true
        count = count + 1;
        try
            Gf = pebi(Gf);
            computeGeometry(Gf);
            break;
        catch
            if count>=4
                flag = 1;
                break;
            end
            h = h/2;
            iter = iter+50;
            if opt.rectangular
                fd = @ (p) drectangle(p, min(xyp(:,1)), max(xyp(:,1)), min(xyp(:,2)), max(xyp(:,2)));
                [p,t] = distmesh_2d(fd, @huniform, h, bbox, iter, ...
                    xyp);
            else
                [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, ...
                    xyp,xyp);
            end
            Gf = triangleGrid(p,t);
        end
    end
end
%--------------------triangle if pebi fails repeatedly--------------------%
if flag == 1
    fprintf('\nCould not successfully generate a PEBI Grid. Switching to triangle grid.\n');
    iter = 50; h = h*2^(count); count = 0; 
    while true
        count = count + 1;
        assert(count<10,' ');
        try
            if opt.rectangular
                fd = @ (p) drectangle(p, min(xyp(:,1)), max(xyp(:,1)), min(xyp(:,2)), max(xyp(:,2)));
                [p,t] = distmesh_2d(fd, @huniform, h, bbox, iter, ...
                    xyp);
            else
                [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, ...
                    xyp,xyp);
            end
            assert(size(t,1)>opt.minTriangles,'');
            close gcf
            Gf = triangleGrid(p,t);
            computeGeometry(Gf);
            computeGeometry(makeLayeredGrid(Gf,1));
            break;
        catch
            h = h/2;
            iter = iter+50;
            
        end
    end
end
%----------------------------gridding complete----------------------------%
nc = [Gf.nodes.coords,repmat(z,Gf.nodes.num,1)];

% xyp = rotatePlane(scaledplane.points,[0 0 1]);
% pdash = [xyp;xyp(1,:)];
% basisCoefficients= bsxfun(@minus,nc,pdash(1,:))*null([0 0 1]);
% nc = bsxfun(@plus,basisCoefficients*null([0 0 1]).', pdash(1,:));

Gf = computeGeometry(makeLayeredGrid(Gf,1));

unitNormal = scaledplane.normal/norm(scaledplane.normal);
tdist = fracplane.aperture/2;

use_points = rotatePlane(nc,+unitNormal);
if any(use_points(:,3)<0)
    use_points = rotatePlane(nc,-unitNormal);
end
dims = max(G.nodes.coords);
use_points = [use_points(:,1)*dims(1) use_points(:,2)*dims(2) use_points(:,3)*dims(3)];
pminus = use_points - repmat(unitNormal,size(use_points,1),1)*tdist;
pplus = use_points + repmat(unitNormal,size(use_points,1),1)*tdist;
nc = [pminus;pplus];
Gf.nodes.coords = nc;
Gf = computeGeometry(Gf);

return

