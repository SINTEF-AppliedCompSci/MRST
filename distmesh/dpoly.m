function d = dpoly ( p, pv )

%*****************************************************************************80
%
%% DPOLY returns the signed distance of one or more points to a polygon.
%
%  Discussion:
%
%    The polygon is described as a sequence of vertices.  In order for
%    a proper calculation, it is necessary that the first vertex be
%    repeated as the last vertex.  Thus, if your polygon is a square,
%    you will be specifying FIVE vertices.
%
%  Licensing:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
%
%  Modified:
%
%    11 March 2006
%
%  Reference:
%
%    Per-Olof Persson and Gilbert Strang,
%    A Simple Mesh Generator in MATLAB,
%    SIAM Review,
%    Volume 46, Number 2, June 2004, pages 329-345.
%
%  Parameters:
%
%    Input, real P(NP,2), the coordinates of a set of nodes.
%
%    Input, real PV(NVS,2), the coordinates of the vertices of the polygon.
%
%    Output, real D, the signed distance of each point to the polygon,
%    which is negative, 0, or positive depending on whether the point
%    is inside, on, or outside the polygon.
%
  np = size ( p, 1 );
  nvs = size ( pv, 1 ) - 1;
%
%  DSEGMENT computes the (unsigned) distance from each point P to every line segment
%  that makes up the polygon.
%
  ds = dsegment ( p, pv );

%ds=zeros(np,nvs);
%for iv=1:nvs
%  ds(:,iv)=donesegment(p,pv(iv:iv+1,:));
%end

  d = min ( ds, [], 2 );
%
%  INPOLYGON is a built-in MATLAB routine, which allows us to determine the
%  sign of the signed distance function.
%
  d = (-1) .^ ( inpolygon ( p(:,1), p(:,2), pv(:,1), pv(:,2) ) ) .* d;

% MEXED

%function ds=donesegment(p,pv)
%
%e=ones(size(p,1),1);
%
%v=diff(pv,1);
%w=p-e*pv(1,:);
%
%c1=sum(w.*v(e,:),2);
%c2=sum(v(e,:).^2,2);
%
%ds=0*e;
%
%ix=c1<=0;
%ds(ix)=sqrt(sum((p(ix,:)-pv(1*ones(sum(ix),1),:)).^2,2));
%
%ix=c1>=c2;
%ds(ix)=sqrt(sum((p(ix,:)-pv(2*ones(sum(ix),1),:)).^2,2));
%
%ix=c1>0 & c2>c1;
%nix=sum(ix);
%if nix>0
%  Pb=ones(nix,1)*pv(1,:)+c1(ix)./c2(ix)*v;
%  ds(ix)=sqrt(sum((p(ix,:)-Pb).^2,2));
%end

  return
end
