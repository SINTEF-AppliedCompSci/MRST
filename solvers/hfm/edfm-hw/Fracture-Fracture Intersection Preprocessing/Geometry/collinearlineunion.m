function line=collinearlineunion(line1,line2)
% pick the furthest points for line1 and line2 and make new line with these
% points. The two lines are assumed to be collinear.

if isempty(line1) && isempty(line2)
    line=[]; return;
elseif isempty(line1) && ~isempty(line2)
    line=line2; return;
elseif ~isempty(line1) && isempty(line2)
    line=line1; return;
end

refline=line1;

refdir=refline(2,:)-refline(1,:);
refdir=refdir/norm(refdir); % make unit
refpoint=refline(1,:);

posline1=line1-repmat(refpoint,2,1);
posline1=dot(posline1,repmat(refdir,2,1),2);

posline2=line2-repmat(refpoint,2,1);
posline2=dot(posline2,repmat(refdir,2,1),2);

posline=[posline1;posline2];
line=[line1;line2];

[~,minind]=min(posline);
[~,maxind]=max(posline);

line=line([minind,maxind],:);

end