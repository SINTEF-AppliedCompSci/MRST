function ds = dsegment ( p, pv )

%*****************************************************************************80
%
%% DSEGMENT computes the distance of points to line segments.
%
%  Discussion:
%
%    This is a pure MATLAB version of the DSEGMENT routine.
%    Normally, the authors of DISTMESH recommend using a faster procedure
%    in which a MATLAB MEX interface is used to call a compiled C code.
%    But if that is not possible, you can fall back on this slower
%    version.
%
%    Note that the line segments are defined by successive pairs of
%    vertices.  In the very common case that the vertices represent
%    a polygon, it is necessary that the list be extended by repeating
%    the first vertex at the end; otherwise the definition of the
%    polygon, and the distance computation, will be incomplete.
%
%  Licensing:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
%
%  Parameters:
%
%    Input, real P(NP,2), a set of points.
%
%    Input, real PV(NVS,2), a set of vertices.  Sequential
%    pairs define NVS-1 line segments.
%
%    Output, real DS(NP,NVS-1); DS(I,J) is the (unsigned) distance
%    of point P(I,:) from line segment PV(J,:)--PV(J+1,:).
%
  nvs = size ( pv, 1 );
  np = size ( p, 1 );

  ds = zeros ( np,nvs-1);

  for iv = 1 : nvs-1
    for ip = 1 : np

      v(1:2) = [ pv(iv+1,1) - pv(iv,1), pv(iv+1,2) - pv(iv,2) ];

      w(1:2) = [ p(ip,1) - pv(iv,1), p(ip,2) - pv(iv,2) ];

      c1 = v * w';
      c2 = v * v';

      if ( c1 <= 0.0 )
        ds(ip,iv) = sqrt ( sum ( ( p(ip,1:2) - pv(iv,1:2) ).^2 ) );
      elseif ( c2 <= c1 )
        ds(ip,iv) = sqrt ( sum ( ( p(ip,1:2) - pv(iv+1,1:2) ).^2 ) );
      else
        ds(ip,iv) = sqrt ( sum ( ( p(ip,1:2) ...
          - pv(iv,1:2) - c1 / c2 * v(1:2) ).^2 ) );
      end

    end 
  end

  return
end
