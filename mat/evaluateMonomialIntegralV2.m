function I = evaluateMonomialIntegralV2(normals, X, Xmid, m)

    I = zeros(1,5);
    Ne = size(normals,1);
    for e = 1:Ne
        nVec = normals(e,:);
        I = I + (1/6*(m(X(e,:)) + m(X(mod(e,Ne)+1,:))) + 2/3*m(Xmid(e,:))).*nVec(1);
    end
    
end