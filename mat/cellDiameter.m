function hK = cellDiameter(X)

    n = size(X,1);
    hK = 0;
    for i = 1:n
        for j = i+1:n
            hK = max(norm(X(i,:)-X(j,:),2),hK);
        end
    end
    
end