function distance = distancePt2Line(p, a, b)
    % p is a matrix of points
    % a is a matrix of starts of line segments
    % b is a matrix of ends of line segments

    pa_distance = sqrt(sum((a-p).^2));
    pb_distance = sqrt(sum((b-p).^2));

    %logical index vector for cases where p is closer to a than b
    idx=pa_distance<pb_distance;

    %assign empty vectors for d and distance
    d=zeros(size(idx));
    distance=zeros(size(idx));

    %compute d for each point
    d(idx)=dot(unit(p(:,idx)-a(:,idx)), unit(b(:,idx)-a(:,idx)));
    d(~idx)=dot(unit(p(:,~idx)-b(:,~idx)), unit(a(:,~idx)-b(:,~idx)));

    %calculate distance
    distance=pb_distance;
    distance(idx)=pa_distance(idx);
    distance(d>0)=distance(d>0).*sqrt(1-d(d>0).^2);
end

function unit_vectors = unit(V)
    % unit computes the unit vectors of a matrix
    % V is the input matrix
    norms = sqrt(sum(V.^2));
    unit_vectors = zeros(size(V));
    normsIndex=norms>eps;
    unit_vectors(:,normsIndex) = V(:,normsIndex)./repmat(norms,size(V,1),1);
end