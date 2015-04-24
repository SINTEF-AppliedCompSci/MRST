function drawAxisCross(s)
if nargin<1, s=.2; end
a = axis;
[xm,ym,zm] = deal(a(1),a(3),a(5));
[xM,yM,zM] = deal(a(2),a(4),a(6));
e = ones(3,1);
v = [s*(xM-xm), s*(yM-ym), s*(zM-zm)];
V = diag(v);
hold on;
quiver3(xm-v(1)*e,ym-v(2)*e,zm-0*e,V(:,1),V(:,2),V(:,3),...
   'LineWidth',2,'Color','k');
text(xm-v(1)*e(1)+V(1,1), ym-v(2)*e(1)+V(1,2), zm+V(1,3), 'x','FontSize',12);
text(xm-v(1)*e(1)+V(2,1), ym-v(2)*e(1)+V(2,2), zm+V(2,3), 'y','FontSize',12);
text(xm-v(1)*e(1)+V(3,1), ym-v(2)*e(1)+V(3,2), zm+V(3,3), 'z','FontSize',12);
hold off
end