function plotPlane(fracplanes, varargin)
if nargin>1
    assert(isnumeric(varargin{1}) & ...
        all(ismember(varargin{1},1:numel(fracplanes))),'Invalid plane index');
    index = varargin{1};
else
    index = 1:numel(fracplanes);
end
for i = 1:numel(index)
    hold on
    points = fracplanes(index(i)).points;
    fill3(points(:,1),points(:,2),points(:,3),rand(1,3));
%     p1 = points(1,:);
%     p2 = points(2,:);
%     p3 = points(3,:);
%     normal = cross(p1-p2, p1-p3);
%     normal = normal./sqrt(sum(normal.^2));
%     A = normal(1); B = normal(2); C = normal(3);
%     D = -dot(normal,p1);
%     xLim = [min(points(:,1)) max(points(:,1))];
%     zLim = [min(points(:,3)) max(points(:,3))];
%     [X,Z] = meshgrid(xLim,zLim);
%     Y = (A * X + C * Z + D)/ (-B);
%     reOrder = [1 2  4 3];
%     patch(X(reOrder),Y(reOrder),Z(reOrder),'b');
%     grid on;
%     alpha(0.3);
end
return
