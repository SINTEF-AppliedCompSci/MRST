function truth=circularfilter(mcellnodes,planenormal,planepoints,tol,varargin)

opt=struct('gaugetol',0,'gaugeradius',0,'center',[],'mcellcentroid',[]);
opt=merge_options(opt,varargin{:});

if isempty(opt.center)
    center=sum(planepoints)/size(planepoints,1);
else
    center=opt.center;
end

if abs(opt.gaugetol)<tol
    gaugetol=norm(max(mcellnodes)-min(mcellnodes));
else
    gaugetol=opt.gaugetol; 
end

if abs(opt.gaugeradius)<tol
    temp=planepoints-repmat(center,size(planepoints,1),1);
    temp=(temp(:,1).^2) + (temp(:,2).^2) +(temp(:,3).^2);
    gaugeradius=realsqrt(max(temp));
else
    gaugeradius=opt.gaugeradius; 
end

if isempty(opt.mcellcentroid)
    mcellcentroid=sum(mcellnodes)/size(mcellnodes,1);
else
    mcellcentroid=opt.mcellcentroid;
end

forbidrad=gaugeradius+gaugetol;

% calculate cell centroid to plane centroid distance along plane.
planenormal=planenormal/norm(planenormal);
relpos=mcellcentroid-center;
normdist=abs(dot(relpos,planenormal));
actdist=norm(relpos);
pardist=realsqrt((actdist^2)-(normdist^2));

if (pardist-forbidrad)>tol % if centroid to center distance is larger than forbidden radius, then no intersection
    truth=false;
else
    truth=true;
end


end