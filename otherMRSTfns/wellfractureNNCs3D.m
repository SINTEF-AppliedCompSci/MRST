function [G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol)
    
    if(isempty(fieldnames(wells))) %Empty wells
        return
    end

    %1. Find the intersection pts between well and fracplanes 
    for wi=1:numel(wells) %number of Wells
        lines=[];
        for si=1:size(wells(wi).points,1)-1 % number of LineSegments
            line=[wells(wi).points(si,:) wells(wi).points(si+1,:)];
            lines=[lines; line];
        end


        XFracPlaneIDs=[];
        XPts=[];
        for pi=1:numel(fracplanes)
            plys=fracplanes(pi).points;
            Pts=LineXPoly3D(lines,plys,tol);
            if ~isempty(Pts)%Found a intersection pts
                XPts=[XPts; Pts];
                XFracPlaneIDs=[XFracPlaneIDs; repmat(pi,size(Pts,1),1)];
            end
        end        
        wells(wi).WellXFracID=XFracPlaneIDs;
        wells(wi).WellXFracPts=XPts;
        fprintf('Total number of fractures that intersect well #%d is %d \n',wi,numel(XFracPlaneIDs));
    end

    %2. Find the fracture cell index
    for wi=1:numel(wells) % number of Wells
        %Find all Mat Cell connects to WellXFracPts
        MatCellIds=Pts2CellIdx(wells(wi).WellXFracPts,G);
        %Find MatCell, FracCell and FracID from MatCellIds, one MatCell  may connect many Fracs
        id=ismember(G.nnc.cells(:,1),MatCellIds);
        MatIdsXFrac=G.nnc.cells(id,1);
        MatXFracCellIds=G.nnc.cells(id,2);
        MatXFracIds=fracCellID2fracID(G,MatXFracCellIds);
        %Find the real WellXFracCellID
        XFracCellIds=zeros(numel(wells(wi).WellXFracID),1);
        for ci=1:numel(wells(wi).WellXFracID)
            idx=find(MatIdsXFrac==MatCellIds(ci));
            idx=idx(MatXFracIds(idx)==wells(wi).WellXFracID(ci));
            XFracCellIds(ci)=MatXFracCellIds(idx);
        end
        %Table of [FracCellID, FracID]
        wells(wi).XFracCellIDs=[XFracCellIds wells(wi).WellXFracID];
    end

    %3. Pass Fracture Aperature into FracGrid
    for i=1:numel(fracplanes)
        G.FracGrid.(['Frac' num2str(i)]).aperture=fracplanes(i).aperture;
    end

end

function xPts = LineXPoly3D(lines,plys,tol)
    %Find the intersection between 3D lines and a 3D polygon
    plane=[plys(1,:),plys(2,:)-plys(1,:),plys(3,:)-plys(1,:)];
    pts=LineXPlane(lines,plane,tol);
    pts = pts(~any(isnan(pts),2),:);
    if ~isempty(pts)
        %Testing if point in the polygon by 2D Projection
        n = size(pts,1);
        Pts_Poly = [pts; plys];
        Pts_Poly2D = ProjectPtsOnPlane(Pts_Poly,plane);
        Pts2D = Pts_Poly2D(1:n,:);
        Poly2D = Pts_Poly2D(n+1:end,:);
        %Matlab 2D inpolygon test
        [in,on] = inpolygon(Pts2D(:,1),Pts2D(:,2),Poly2D(:,1),Poly2D(:,2));
        xPts=pts(in|on,:);
    else
        xPts=[];
    end
    
end

function pts = LineXPlane(lns,pln,tol)
    % intersection between lines and a plane, 3D                                    
    n = size(lns,1);
    nrm = repmat(cross(pln(4:6),pln(7:9),2),n,1);
    lns = [lns(:,1:3),lns(:,4:6)-lns(:,1:3)];
    f = (abs(dot(nrm,lns(:,4:6),2)) >= tol);
    pts = nan(n,3);
    d = bsxfun(@minus,pln(1:3),lns(:,1:3));
    t = dot(nrm(f,:),d(f,:),2)./dot(nrm(f,:),lns(f,4:6),2);
    r = (-tol < t) & (t < (1+tol));
    f(f) = r;
    %pts(f,:) = lns(f,1:3)+t(r).*lns(f,4:6);
    if(~isempty(t(r)))
        pts(f,:) = lns(f,1:3)+bsxfun(@times,t(r),lns(f,4:6));
    end
end

function [ots] = ProjectPtsOnPlane(pts,pln)
    %project 3D points on target plane
    vec = @(v)sqrt(sum(v.*v,2));
    
    p0 = pln(1:3);
    d1 = pln(4:6);
    d2 = pln(7:9);
    n = size(pts,1);
    if n == 1                                                                   % only one point
        s = dot(bsxfun(@minus,pts,p0),d1,2)./vec(d1);
        t = dot(bsxfun(@minus,pts,p0),d2,2)./vec(d2);
    else                                                                        % many points
        d1 = d1/vec(d1);
        d2 = d2/vec(d2);
        f = ones(n,1);
        s = dot(bsxfun(@minus,pts,p0),d1(f,:),2);
        t = dot(bsxfun(@minus,pts,p0),d2(f,:),2);
    end
    ots = [s,t];
end

function FracIds=fracCellID2fracID(G,FracCellIds)
    %Convert fracture cell ID into fracplane IDs
    NumFracs=length(fieldnames(G.FracGrid));
    
    FracIds=zeros(numel(FracCellIds),1);
    for fi=1:NumFracs
        startID=G.FracGrid.(['Frac' num2str(fi)]).cells.start;
        endID=startID+G.FracGrid.(['Frac' num2str(fi)]).cells.num;
        FracIds(FracCellIds>= startID & FracCellIds< endID)=fi;
    end
    
end