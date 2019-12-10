function cub = compressCubature(cub, type, deg, compressionType)
    % Compress the cubature cub defined on the entitiy_type
    % ('volume' or 'face'). The options are named to match
    % getCubature.
    %
    % The cubature is valid to degree deg, i.e. integrals of this
    % polynomial order is integrated exactly.
    %
    % compression_type is either 'Tchakaloff' (uses nonlinear least
    % squares) or 'Fekete' (uses qr factorization with column
    % pivoting). Default is 'Fekete'.

    % http://www.math.unipd.it/~marcov/mysoft/subp/comprexcub.m
    % Compression of bivariate cubature formulas by Tchakaloff points
    % or approximate Fekete points
    % useful, for example, in node reduction of algebraic cubature formulas
    % see the web page: http://www.math.unipd.it/~marcov/Tchakaloff.html

    % by Federico Piazzon, Alvise Sommariva and Marco Vianello
    % University of Padova, May 2016

    if nargin == 3
        compressionType = 'Fekete';
    end

    assert(strcmp(compressionType, 'Fekete') |...
           strcmp(compressionType, 'Tchakaloff'))

    G = cub.G;
    switch type
        case 'volume'
            num_entities = G.cells.num;
        case 'face'
            num_entities = G.faces.num;
        otherwise
            error(['Unknown entities ', type])
    end

    % The compressed rule contains ncub points, where ncub is the
    % dimension of the polynomial space. 
    ncompcub = nchoosek(deg+cub.dim, cub.dim);

    % Only compress the elements with many points (for now we only
    % assume it's all or none)
    entities = 1:num_entities;
    dp = diff(cub.pos);
    entities = entities(dp > ncompcub);
    assert(numel(entities) == num_entities | ...
        numel(entities) == 0);
    if numel(entities) == 0
        return;
    end
    
    points  = zeros(G.cells.num*ncompcub, G.griddim);
    weights = zeros(G.cells.num*ncompcub, 1);

    for e = entities
        % Get the cubature rule for this cell
        [W, X, omega] = getCubature(cub, e, type);

        % total-degree Chebyshev-Vandermonde matrix at X
        rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
        V=chebvand(deg,X,rect);
        % polynomial basis orthogonalization
        [Q,R]=qr(V,0);
        % tiny imaginary parts could appear
        Q=real(Q);
        % possible re-orthogonalization
        % [Q,R]=qr(Q,0);

        % moments of the orthogonal basis
        orthmom=Q'*omega;
        % weigths computation
        switch compressionType
            case 'Tchakaloff'
                % Tchakaloff points (positive weights)
                weights = lsqnonneg(Q',orthmom);
                %[weights, resnorm, residual]=lsqnonneg(Q',orthmom);
                % C=Q';
                % d=orthmom;
                % 10*max(size(C))*norm(C,1)*eps
                % resnorm
                % residual
            case 'Fekete'
                % approximate Fekete points (possible negative weights)
                weights=Q'\orthmom;
        end
        % indexes of nonvanishing weights and compression
        ind=find(abs(weights)>0);
        %pts=X(ind,:);
        %w=weights(ind);
        %[size(omega), size(w)]
        % moment reconstruction error
        % bivariate Chebyshev basis
        % mom=V'*omega;
        % momerr=norm(V(ind,:)'*w-mom);
        % discrete OP basis
        %momerr=norm(Q(ind,:)'*w-orthmom);

        % Update cub
        % NB: We assume we compress all entities
        idx = (ncompcub*(e-1)+1):ncompcub*e;
        points(idx, :) = X(ind, :);
        weights(idx) = weights(ind);
    end

    % Update cub data
    % NB: We assume we compress all entities
    cub.pos(2:end) = cumsum(ncompcub*ones(num_entities, 1)) + 1;
    cub.points = points;
    cub.weights = weights;
    
end

function V = chebvand(deg,X,rect)

    % INPUT:
    % deg = polynomial degree
    % X = 2-column array of point coordinates
    % rect = 4-component vector such that the rectangle
    % [rect(1),rect(2)] x [rect(3),rect(4)] contains X

    % OUTPUT:
    % V = Chebyshev-Vandermonde matrix at X, graded lexic. order

    % FUNCTION BODY
    % rectangle containing the mesh
    if isempty(rect)
        rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
    end

    % couples with length less or equal to deg
    % graded lexicographical order
    j=(0:1:deg);
    [j1,j2]=meshgrid(j);
    dim=(deg+1)*(deg+2)/2;
    couples=zeros(dim,2);
    for s=0:deg
        good=find(j1(:)+j2(:)==s);
        couples(1+s*(s+1)/2:(s+1)*(s+2)/2,:)=[j1(good) j2(good)];
    end

    % mapping the mesh in the square [-1,1]^2
    a=rect(1);b=rect(2);c=rect(3);d=rect(4);
    map=[(2*X(:,1)-b-a)/(b-a) (2*X(:,2)-d-c)/(d-c)];

    % Chebyshev-Vandermonde matrix on the mesh
    T1=chebpolys(deg,map(:,1));
    T2=chebpolys(deg,map(:,2));
    V=T1(:,couples(:,1)+1).*T2(:,couples(:,2)+1);

end



function T=chebpolys(deg,x)

    % computes the Chebyshev-Vandermonde matrix on the real line by recurrence

    % INPUT:
    % deg = maximum polynomial degree
    % x = 1-column array of abscissas

    % OUTPUT
    % T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg

    T=zeros(length(x),deg+1);
    t0=ones(length(x),1);
    T(:,1)=t0;
    t1=x;
    T(:,2)=t1;

    for j=2:deg
        t2=2*x.*t1-t0;
        T(:,j+1)=t2;
        t0=t1;
        t1=t2;
    end

end