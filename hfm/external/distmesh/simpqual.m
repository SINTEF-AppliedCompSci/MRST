function q = simpqual ( p, t, type )

%*****************************************************************************80
%
%% SIMPQUAL computes the simplex quality.
%
%  Discussion:
%
%    If TYPE = 1, the quality measure used computes the ratio of
%    the radii of the inscribed and circumscribed circles or spheres.
%
%    Only dimensions 1, 2 and 3 can be handled.  In 1D, the quality is always 1.
%
%  Licensing:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
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
%    Input, real P(NP,ND), the coordinates of a set of nodes.
%
%    Input, integer T(NT,1:ND+1), a list of the nodes which make up each 
%    simplex in the mesh.
%
%    Input, integer TYPE, specifies the quality measure:
%    1: Compute the radius ratio (default)
%    2: Use an approximate formula.
%
%    Output, real Q(NT), the simplex quality measure.
%
  if ( nargin < 3 )
    type = 1;
  end

  switch type
%
%  RADIUS RATIO
%
    case 1

      switch size ( p, 2 )

        case 1

          q = ones ( 1, size(t,2) );

        case 2

          a=sqrt(sum((p(t(:,2),:)-p(t(:,1),:)).^2,2));
          b=sqrt(sum((p(t(:,3),:)-p(t(:,1),:)).^2,2));
          c=sqrt(sum((p(t(:,3),:)-p(t(:,2),:)).^2,2));
          r=1/2*sqrt((b+c-a).*(c+a-b).*(a+b-c)./(a+b+c));
          R=a.*b.*c./sqrt((a+b+c).*(b+c-a).*(c+a-b).*(a+b-c));
          q = 2 * r ./ R;

        case 3

          d12=p(t(:,2),:)-p(t(:,1),:);
          d13=p(t(:,3),:)-p(t(:,1),:);
          d14=p(t(:,4),:)-p(t(:,1),:);
          d23=p(t(:,3),:)-p(t(:,2),:);
          d24=p(t(:,4),:)-p(t(:,2),:);
          d34=p(t(:,4),:)-p(t(:,3),:);
          v=abs(dot(cross(d12,d13,2),d14,2))/6;
          s1=sqrt(sum(cross(d12,d13,2).^2,2))/2;
          s2=sqrt(sum(cross(d12,d14,2).^2,2))/2;
          s3=sqrt(sum(cross(d13,d14,2).^2,2))/2;
          s4=sqrt(sum(cross(d23,d24,2).^2,2))/2;
          p1=sqrt(sum(d12.^2,2)).*sqrt(sum(d34.^2,2));
          p2=sqrt(sum(d23.^2,2)).*sqrt(sum(d14.^2,2));
          p3=sqrt(sum(d13.^2,2)).*sqrt(sum(d24.^2,2));
          q=216*v.^2./(s1+s2+s3+s4)./sqrt((p1+p2+p3).*(p1+p2-p3).* ...
                                    (p1+p3-p2).*(p2+p3-p1));

        otherwise

          error ( 'SIMPQUAL - Dimension not implemented.' );
    end
%
%  APPROXIMATE FORMULA
%
    case 2

     switch size(p,2)

       case 1

         q = ones ( 1, size(t,2) );

       case 2

         d12 = sum((p(t(:,2),:)-p(t(:,1),:)).^2,2);
         d13 = sum((p(t(:,3),:)-p(t(:,1),:)).^2,2);
         d23 = sum((p(t(:,3),:)-p(t(:,2),:)).^2,2);
         q = 4*sqrt(3)*abs(simpvol(p,t))./(d12+d13+d23);

       case 3

         d12 = sum((p(t(:,2),:)-p(t(:,1),:)).^2,2);
         d13 = sum((p(t(:,3),:)-p(t(:,1),:)).^2,2);
         d14 = sum((p(t(:,4),:)-p(t(:,1),:)).^2,2);
         d23 = sum((p(t(:,3),:)-p(t(:,2),:)).^2,2);
         d24 = sum((p(t(:,4),:)-p(t(:,2),:)).^2,2);
         d34 = sum((p(t(:,4),:)-p(t(:,3),:)).^2,2);
         q = 216*abs(simpvol(p,t))/sqrt(3)./(d12+d13+d14+d23+d24+d34).^(3/2);

       otherwise

         error ( 'SIMPQUAL - Dimension not implemented.' );
    end

  otherwise

    error ( 'SIMPQUAL - Incorrect type.' );

  end

  return
end
