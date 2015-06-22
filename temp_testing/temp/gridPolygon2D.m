function nodes = gridPolygon2D(xyp,varargin)

% opt = struct('ratio', 0.5);
% opt = merge_options(opt, varargin{:});

possible_ratios = 0:0.0001:1;
possible_ratios = possible_ratios(mod(1,possible_ratios)==0);
%
div = 3;
pfull = [xyp;xyp(1,:)];
[h,side] = min(sqrt(sum(diff(pfull,1).^2,2))); % min side length
hmax = max(sqrt(sum(diff(pfull,1).^2,2)));
ratio = h/div;
[~,loc] = min(abs(possible_ratios-ratio));
ratio = possible_ratios(loc);

t = transpose(0:ratio:1);
x = pfull(side,1)*(1-t) + pfull(side+1,1)*t;
y = pfull(side,2)*(1-t) + pfull(side+1,2)*t;
gridp = [x,y];
minl = h*div/(div-1);
% minr = zeros(size(xyp,1),1);
for i = 1:size(xyp,1)
    pts = pfull(i:i+1,:);
    if i~=side
        l = sqrt(sum(diff(pts,1).^2));
        r = minl/l;
        [~,loc] = min(abs(possible_ratios-r));
        r = possible_ratios(loc);
        t = transpose(0:r:1);
        x = pts(1,1)*(1-t) + pts(2,1)*t;
        y = pts(1,2)*(1-t) + pts(2,2)*t;
        gridp = [gridp;x,y]; %#ok
%         minr(i) = r;
    end
end
% move = minr;
% meanp = mean(xyp,1);
xypshift(:,1) = xyp(:,1);%*(1-move) + move*(meanp(:,1));
xypshift(:,2) = xyp(:,2);%*(1-move) + move*(meanp(:,2));
minxy = min(xypshift);
maxxy = max(xypshift);
r = minl/sqrt(sum((maxxy-minxy).^2));
[X,Y] = meshgrid(linspace(minxy(1),maxxy(1),1/r),linspace(minxy(2),maxxy(2),1/r));
padd = [X(:) Y(:)];
padd = padd(inpolygon(padd(:,1),padd(:,2),xypshift(:,1),xypshift(:,2)),:);
ldist = zeros(size(padd,1),size(xyp,1));
diffxy = diff(pfull,1);
m = diffxy(:,2)./diffxy(:,1);
for i = 1:size(xyp,1)
    c = -m(i)*xyp(i,1) + xyp(i,2);
    ldist(:,i) = abs((m(i)*padd(:,1) - padd(:,2) + c)./sqrt(1 + m(i)^2));
end
distdiff = ldist>=minl/5;
padd = padd(sum(distdiff,2)==size(xyp,1),:);
nodes = unique([xyp;gridp;padd],'rows','stable');
return
        

