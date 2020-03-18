function [A, b, I, J, V] = AssemAb(G, ~, mu, rho, T, W, src, well_posed)
dispif(mrstVerbose, 'AssemAb\n');
ncf = max(diff(G.cells.facePos));
nc = G.cells.num;

num_W_rate = getNoRateWells(W);
nc = nc + num_W_rate;
b = 0 * repmat(T{1}(1), nc, 1);
k = 1;
%[I,J,V]=deal(zeros(ncf*nc,1));k=1;
[I, J] = deal(zeros(ncf*nc, 1));
V = 0 * repmat(T{1}(1), ncf*nc, 1);

for i_face = 1:G.faces.num
    c1 = G.faces.neighbors(i_face, 1);
    c2 = G.faces.neighbors(i_face, 2);
    if (all([c1, c2] ~= 0))
        I(k) = c1;
        J(k) = c1;
        V(k) = T{1}(i_face);
        k = k + 1;
        I(k) = c1;
        J(k) = c2;
        V(k) = -T{2}(i_face);
        k = k + 1;
        I(k) = c2;
        J(k) = c2;
        V(k) = T{2}(i_face);
        k = k + 1;
        I(k) = c2;
        J(k) = c1;
        V(k) = -T{1}(i_face);
        k = k + 1;
        % I(k)=c1;J(k)=c1;V(k)=value(T{1}(i_face));k=k+1;
        % I(k)=c1;J(k)=c2;V(k)=value(-T{2}(i_face));k=k+1;
        % I(k)=c2;J(k)=c2;V(k)=value(T{2}(i_face));k=k+1;
        % I(k)=c2;J(k)=c1;V(k)=value(-T{1}(i_face));k=k+1;
    else
        c1 = max(c1, c2);
        I(k) = c1;
        J(k) = c1;
        V(k) = T{1}(i_face);
        %V(k)=value(T{1}(i_face));
        k = k + 1;
        %T2_val = T{2}(i_face) + 0*V(1);
        %b(c1) = b(c1) + T2_val;
        b(c1) = b(c1) + T{2}(i_face);
    end
end

%----------------------------------------------------------
g = gravity();
for i = 1:numel(W)
    if (strcmpi(W(i).type, 'bhp'))
        pbh = W(i).val;
        dZ = W(i).dZ;
        for j = 1:numel(W(i).cells)
            mycell = W(i).cells(j);
            I(k) = mycell;
            J(k) = mycell;
            V(k) = W(i).WI(j) / mu;
            k = k + 1;
            b(mycell) = b(mycell) + W(i).WI(j) / mu * (pbh + rho * g(3) * dZ(j));
        end
    elseif (strcmpi(W(i).type, 'rate'))
        rate = W(i).val;
        dZ = W(i).dZ;
        for j = 1:numel(W(i).cells)
            mycell = W(i).cells(j);
            I(k) = mycell;
            J(k) = mycell;
            V(k) = W(i).WI(j) / mu;
            k = k + 1;
            I(k) = mycell;
            J(k) = G.cells.num + i;
            V(k) = -W(i).WI(j) / mu;
            k = k + 1;
            I(k) = G.cells.num + i;
            J(k) = G.cells.num + i;
            V(k) = W(i).WI(j) / mu;
            k = k + 1;
            I(k) = G.cells.num + i;
            J(k) = mycell;
            V(k) = -W(i).WI(j) / mu;
            k = k + 1;
            b(mycell) = b(mycell) + W(i).WI(j) / mu * (rho * g(3) * dZ(j));
            b(G.cells.num+i) = b(G.cells.num+i) + W(i).WI(j) / mu * (rho * g(3) * dZ(j));
        end
        b(G.cells.num+i) = b(G.cells.num+i) + rate;
    else
        % write code here bababababababbababaabababababababababababababba
        error('code under development!')
    end
end
%-------------------------------------------------------
I(k:end) = [];
J(k:end) = [];
V(k:end) = [];

if isa(V, 'ADI')
    A = [];
    b = [];
    return;
else
    A = sparse(I, J, V, nc, nc);
end

if ~isempty(src)
    b(src.cell(:)) = b(src.cell(:)) + src.rate(:);
end

% Fix well-posedness similar to incompTPFA
if ~well_posed
    if A(1) > 0
        A(1) = 2 * A(1);
    else
        keyboard
    end
end

%keyboard
end
