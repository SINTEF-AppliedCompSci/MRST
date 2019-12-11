function clippoints=polygonplaneclip(points,planepoint,planedirection,tol)
% 3D polygon vs plane clipping tool. This function takes in the points of a
% polygon (arranged!) in the form of an Nx3 matrix, each row containing xyz
% coordinates of the vertices. The plane is given by a point on the plane
% (row vector) and a plane direction. Take note that polygon points on the
% positive side of the plane are retained. tol is a floating number
% tolerance.

if isempty(points)  % if input points is null set. Return null set.
    clippoints=points;
    return;
end

planedirection=planedirection/norm(planedirection); % make unit normal

numpoints=size(points,1);
relpos=points-repmat(planepoint,numpoints,1);

normaldist=dot(relpos,repmat(planedirection,numpoints,1),2);

sign=normaldist; 
sign(sign>tol)=1; 
sign(sign<-tol)=-1; 
sign(abs(sign)<tol)=0; % keep points that are on plane
sign=[sign;sign(1,:)];
points=[points;points(1,:)];

clippoints=-1*ones(numpoints,3); % preallocate upper bound, remove excess later
count=0;

for i=1:numpoints
    if sign(i)==1 || sign(i)==0 % save point i if it is on positive side or on plane
        count=count+1;
        clippoints(count,:)=points(i,:);
    end
    
    diff=sign(i+1)-sign(i); % sign change (negative means next step crosses into negative zone, vice versa)
    % evaluate next point
    if diff==0 % no sign change, no intersection to worry about.
        continue;
    else % sign change (going into a new zone)
        if sign(i)==0
            continue; % currently on plane, no worries about intersection, next loop will save this point
        elseif sign(i)==1 && sign(i+1)~=0 % next point crosses into negative zone from positive
            % perform intersection calc
            linepts=points(i:(i+1),:);
            [~,xpoint,~] = linesegmentplaneintersect(linepts,planepoint,planedirection,tol);
            % replace clippoint(i+1) with intersection since next point would be
            % ignored in next iteration (sign(i+1)==-1)
            count=count+1;
            clippoints(count,:)=xpoint;
        elseif sign(i)==-1 && sign(i+1)~=0 % next point crosses into positive zone from negative
            % perform intersection calc
            linepts=points(i:(i+1),:);
            [~,xpoint,~] = linesegmentplaneintersect(linepts,planepoint,planedirection,tol);
            % replace clippoint(i) with intersection since the current
            % point(i) has been ignored
            count=count+1;
            clippoints(count,:)=xpoint;
        elseif sign(i+1)==0
            % if there is a sign change and next point is on plane. Then
            % don't do anything. The next for loop will save the point on
            % plane
            continue;
        end
    end
end

clippoints=removeexcess(clippoints,-1);

% final output should be arranged sequentially

end