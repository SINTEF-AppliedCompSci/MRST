function hK = cellDiameter(X)

    n = size(X,1);
    hK = 0;
    for i = 1:n
        hK = max(norm(repmat(X(i,:),n,1)-X),hK);
    end
    
end