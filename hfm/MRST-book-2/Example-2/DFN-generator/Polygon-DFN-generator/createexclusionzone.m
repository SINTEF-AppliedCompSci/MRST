function fracture=createexclusionzone(fracture,allowance,varargin)
% given points that describe a convex polygon, this function defines an
% exclusion zone using the allowance. The exclusion zone is an extruded
% expanded polygon. The polygon vertices are moved further away by
% allowance*distance and the thickness of the exclusion zone is
% allowance*meandistance well.
%
% the output will be a struct class object containing the fields surfaces.
% Each field will contain information on the vertices of the surfaces and
% the normal direction of the surfaces pointing into the exclusion zone.
% The surfaces together form the exclusion zone.
%
% An option is provided to triangulate all surfaces for later intersection
% checking. Triangulation will create a triangle list for every surface.
% Each row in this triangle list will contain indices of nodes that form a
% triangle.

opt=struct('triangulate',true,'visualize',false,'arranged',true);
opt=merge_options(opt,varargin{:});

%% Create exclusion zone
fsize=fracture.size;
points=fracture.points;
numpoints=size(points,1);

planenormal=fracture.normal;

if ~opt.arranged
    points=arrangenodes(points,planenormal);
end

centroid=sum(points)/numpoints;

tempdirs=points-repmat(centroid,numpoints,1);
tempdirslength=realsqrt((tempdirs(:,1).^2)+(tempdirs(:,2).^2)+(tempdirs(:,3).^2));
tempdirs=tempdirs./repmat(tempdirslength,1,3);
extendpoints=points+allowance*fsize*tempdirs;

% surface 1 calculation: one side of the polygon
% meandist=points-repmat(centroid,numpoints,1);
% meandist=realsqrt((meandist(:,1).^2)+(meandist(:,2).^2)+(meandist(:,3).^2));
% meandist=sum(meandist)/numpoints;
exclusionzone.surface1.points=extendpoints+0.5*allowance*fsize*repmat(planenormal,numpoints,1);
exclusionzone.surface1.insidedirection=-planenormal;

% surface 2 calculation: other side of the polygon
exclusionzone.surface2.points=extendpoints-0.5*allowance*fsize*repmat(planenormal,numpoints,1);
exclusionzone.surface2.insidedirection=planenormal;

% here on out, generate the side surfaces of the exclusion zone. Number of
% surfaces = number of points
surfpoints1=[exclusionzone.surface1.points;exclusionzone.surface1.points(1,:)];
surfpoints2=[exclusionzone.surface2.points;exclusionzone.surface2.points(1,:)];
for i=1:numpoints
    fieldname=['surface',num2str(i+2)];
    % use surface 1 and surface 2 points to define side surfaces.
    % Arrangement is guaranteed through the way the points are saved.
    surfpointsi=[surfpoints1(i:(i+1),:);surfpoints2([(i+1),i],:)];
    normali=cross(planenormal,surfpointsi(2,:)-surfpointsi(4,:));
    normali=normali/norm(normali);
    if dot(centroid-surfpointsi(1,:),normali)<0 
        % making sure that the normal points inwards
        normali=-normali;
    end
    exclusionzone.(fieldname).points=surfpointsi;
    exclusionzone.(fieldname).insidedirection=normali;    
end

%% Triangulation of surfaces
if opt.triangulate % triangulate if required by user
    numsurfaces=numpoints+2;
    for i=1:numsurfaces
        fieldname=['surface',num2str(i)];
        surfpointsi=exclusionzone.(fieldname).points;
        numpointsi=size(surfpointsi,1);
        numtrianglesi=numpointsi-2;
        % points are already arranged.
        trilisti=[ones(numtrianglesi,1),(2:(numpointsi-1))',(3:numpointsi)'];
        exclusionzone.(fieldname).TriangleList=trilisti;
    end
end

%% Plot if required
if opt.visualize && ~opt.triangulate
    figure;
    numsurfaces=numpoints+2;
    for i=1:numsurfaces
        fieldname=['surface',num2str(i)];
        X=exclusionzone.(fieldname).points(:,1);
        Y=exclusionzone.(fieldname).points(:,2);
        Z=exclusionzone.(fieldname).points(:,3);
        fill3(X,Y,Z,'y');
        alpha(0.5);
        hold on;
    end
    X=points(:,1);
    Y=points(:,2);
    Z=points(:,3);
    fill3(X,Y,Z,'b');
    set(gca,'Zdir','reverse')
    view(30,30);
    axis equal tight;
    hold off;
elseif opt.visualize && opt.triangulate
    figure; hold on;
    numsurfaces=numpoints+2;
    for i=1:numsurfaces
        fieldname=['surface',num2str(i)];
        trilist=exclusionzone.(fieldname).TriangleList;
        X=exclusionzone.(fieldname).points(:,1);
        Y=exclusionzone.(fieldname).points(:,2);
        Z=exclusionzone.(fieldname).points(:,3);
        for j=1:size(trilist,1);
            fill3(X(trilist(j,:)),Y(trilist(j,:)),Z(trilist(j,:)),'y');
        end
        alpha(0.5);
    end
    X=points(:,1);
    Y=points(:,2);
    Z=points(:,3);
    fill3(X,Y,Z,'b');
    set(gca,'Zdir','reverse')
    view(30,30);
    axis equal tight;
    hold off;
end

fracture.exclusionzone=exclusionzone;
end