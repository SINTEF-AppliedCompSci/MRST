function T = TransNTPFA(G, bc, mu, u, OSflux)
    dispif(mrstVerbose, 'TransNTPFA\n');

    if numel(u) > G.cells.num
        % FlowNTPFA compability
        u = u(1:G.cells.num);
    end

    m = G.faces.num;
    n = G.cells.num;
    t1end = zeros(m, 1);
    t2end = zeros(m, 1);
    t11 = zeros(m, 2);
    t22 = zeros(m, 2);

    nij1 = zeros(m, 1);
    nij2 = zeros(m, 1);

    internal = (1:m);
    internal(~all(G.faces.neighbors ~= 0, 2)) = [];

    % Sweep for sparsity set up
    for i = internal
        t1 = OSflux{i, 1};
        t2 = OSflux{i, 2};
        nij1(i) = max(0, size(t1, 1)-3);
        nij2(i) = max(0, size(t2, 1)-3);
        t11(i, 1:2) = [t1(1, 2), t1(2, 2)];
        t22(i, 1:2) = [t2(1, 2), t2(2, 2)];
    end

    % Set up sparsity pattern and values
    ii1 = zeros(sum(nij1), 1);
    ii2 = zeros(sum(nij2), 1);
    jj1 = ii1;
    jj2 = ii2;
    vv1 = ii1;
    vv2 = ii2;
    s1 = [1; cumsum(nij1) + 1];
    s2 = [1; cumsum(nij2) + 1];

    % Don't loop over zero rows
    nzrows1 = (1:m);
    nzrows2 = (1:m);
    nzrows1(nij1 == 0) = [];
    nzrows2(nij2 == 0) = [];

    for i = nzrows1
        t1 = OSflux{i, 1};
        idx1 = s1(i):(s1(i+1) - 1);
        ii1(idx1) = i;
        jj1(idx1) = t1(3:end-1, 1);
        vv1(idx1) = t1(3:end-1, 2);
        t1end(i) = t1(end, 2);
    end
    
    for i = nzrows2
        t2 = OSflux{i, 2};
        idx2 = s2(i):(s2(i+1) - 1);
        ii2(idx2) = i;
        jj2(idx2) = t2(3:end-1, 1);
        vv2(idx2) = t2(3:end-1, 2);
        t2end(i) = t2(end, 2);
    end

    t1sp = sparse(ii1, jj1, vv1, m, n);
    t2sp = sparse(ii2, jj2, vv2, m, n);

    r1 = t1sp * u + t1end;
    r2 = t2sp * u + t2end;

    epstol = 1e-12 * max(full(max(t1sp, [], 2)), full(max(t2sp, [], 2)));
    ir1 = abs(r1) <= epstol;
    ir2 = abs(r2) <= epstol;
    r1(ir1) = 0;
    r2(ir2) = 0;

    jj = abs(r1+r2) > epstol;
    mu1 = 0 * r1 + 0.5;
    mu2 = 0 * r2 + 0.5;
    mu1(jj) = r2(jj) ./ (r1(jj) + r2(jj));
    mu2(jj) = ones(sum(jj), 1) - mu1(jj);

    T = cell(2, 1);
    T{1} = 0 * r1;
    T{2} = 0 * r2;
    T{1}(internal) = mu1(internal) .* t11(internal, 1) + mu2(internal) .* t22(internal, 2);
    T{2}(internal) = mu1(internal) .* t11(internal, 2) + mu2(internal) .* t22(internal, 1);

    %     TT = T;
    %
    %     T = cell(2, 1);
    %     T{1} = 0 * repmat(u(1), G.faces.num, 1);
    %     T{2} = 0 * repmat(u(1), G.faces.num, 1);
    %     for i_face = 1:G.faces.num
    %         if all(G.faces.neighbors(i_face, :) ~= 0)
    %             t1 = OSflux{i_face, 1};
    %             t2 = OSflux{i_face, 2};
    %             r1 = t1(3:end-1, 2)' * u(t1(3:end-1, 1)) + t1(end, 2);
    %             r2 = t2(3:end-1, 2)' * u(t2(3:end-1, 1)) + t2(end, 2);
    %             eps = 0 * r1 + 1e-12 * max(abs([t1(:, end); t2(:, end)]));
    %             if (abs(r1) <= eps)
    %                 r1 = 0 * r1;
    %             end
    %             if (abs(r2) <= eps)
    %                 r2 = 0 * r2;
    %             end
    %
    %             if (abs(r1+r2) > eps)
    %                 mu1 = r2 ./ (r1 + r2); % ./ for AD
    %                 mu2 = 1 - mu1;
    %             else
    %                 mu1 = 0 * r1 + 0.5;
    %                 mu2 = 0 * r2 + 0.5;
    %             end
    %             T{1}(i_face) = (mu1 * t1(1, 2) + mu2 * t2(2, 2));% ./ mu(i_face);
    %             T{2}(i_face) = (mu1 * t1(2, 2) + mu2 * t2(1, 2));% ./ mu(i_face);
    %         else
    % %             ind = find(bc.face == i_face, 1);
    % %             if strcmpi(bc.type{ind}, 'pressure')
    % %                 t1 = OSflux{i_face, 1};
    % %                 t2 = OSflux{i_face, 2};
    % %                 t11 = t1(1, 2);
    % %                 t12 = t1(2, 2);
    % %                 t22 = t2(1, 2);
    % %                 t21 = t2(2, 2);
    % %                 r1 = t1(3:end-1, 2)' * u(t1(3:end-1, 1)) + t1(end, 2);
    % %                 r2 = t2(end, 2) + 0 .* u(1);
    % %                 eps = 0 * r1 + 1e-12 * max(abs([t1(:, end); t2(:, end)]));
    % %                 if (abs(r1) <= eps), r1 = 0; end
    % %                 if (abs(r2) <= eps), r2 = 0; end
    % %                 if (abs(r1+r2) > eps)
    % %                     mu1 = r2 ./ (r1 + r2);
    % %                     mu2 = 1 - mu1;
    % %                 else
    % %                     mu1 = 0.5;
    % %                     mu2 = 0.5;
    % %                 end
    % %                 T{1}(i_face) = mu1 * t11 + mu2 * t21; % Divide by mu as above?
    % %                 T{2}(i_face) = (mu1 * t12 + mu2 * t22) * bc.value{ind}(G.faces.centroids(i_face, :));
    % %             else
    % %                 T{2}(i_face) = bc.value{ind}(G.faces.centroids(i_face, :));
    % %             end
    %         end
    %     end
    %
    %     % check
    %     for i = 1:G.faces.num
    %         if all(G.faces.neighbors(i, :) ~= 0)
    %             e1 = abs(T{1}(i)-TT{1}(i));
    %             e2 = abs(T{2}(i)-TT{2}(i));
    %             if e1 > 0 || e2 > 0
    %                 T{1}(i), TT{1}(i), e1
    %                 T{2}(i), TT{2}(i), e2
    %                 keyboard
    %             end
    %         end
    %     end
end
