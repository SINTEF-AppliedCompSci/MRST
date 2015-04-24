function [xOpt, poly] = argmaxCubic(p1, p2)
% find max of cubic polynomial through p1, p2
% shift values:
shift = p1.a;
p1.a = 0;
p2.a = p2.a-shift;
poly    = zeros(4,1);
poly(3:4) = [p1.dv; p1.v];
a = p2.a;
poly(1:2) = [a^3, a^2; 3*a^2 2*a]\[p2.v-p1.dv*a-p1.v; p2.dv-p1.dv];
xe = roots(polyder(poly));
if numel(xe) == 0
    xe = -inf;
elseif any(imag(xe)~=0) 
    xe = inf;
elseif xe(1) == xe(end) % single root
    if poly(1) ~= 0
        xe = inf;
    elseif poly(2) < 0
        xe = xe(1);
    else
        xe = -inf;
    end
else  % two roots
    [~,ix] = max(polyval(poly,xe));
    xe = xe(ix);
    if xe < p1.a
        xe = -inf;
    end
end
xOpt = xe + shift;
    
%--------------------------------------------------------------------------
figure(99), clf, hold on
x1 = p1.a; x2 = p2.a;
xi = x1/2:(x2-x1)/50:2*x2;
plot(xi+shift, polyval(poly, xi));
plot(xe+shift, polyval(poly, xe), 'o', 'Markersize', 14, 'MarkerEdgeColor', 'r');
plot(p1.a+shift, p1.v, 'o', 'Markersize', 14);
plot(p2.a+shift, p2.v, 'o', 'Markersize', 14);
l0 = p1.dv*xi + p1.v -p1.dv*p1.a;
lx = p2.dv*xi + p2.v -p2.dv*p2.a;
plot(xi+shift, l0, '--k')
plot(xi+shift, lx, '--k')
drawnow
end

