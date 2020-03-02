function out = MPSA_vectorized(G,C,varargin)

% TODO
% Proper handling of Neumann Boundaries
% Stabilization term

opt=struct('method','O',...
    'bc',[],...
    'weakCont',1,'eta',0,'biots',0);

opt=merge_options(opt,varargin{:});
bc=opt.bc;
if ~strcmpi(opt.method,'O')
    error('Have not implemented other options than O-method')
end

% Bookkeeping
Nd = G.griddim;
Nc = G.cells.num;
Nf = G.faces.num;
Nn = G.nodes.num;

% This is meant to introduce a notion of active cells, can be of use in
% plasticity etc. Implementation should be straightforward (find active
% nodes, local solves only for these). Note that combining old and new
% half-stressmissibilities may need some care.
% if nargin == 2 || numel(ci) == 0
%     ci = 1 : Nc;
% else
%     ci = reshape(ci,1,[]);
% end

%% Identifiers for boundary conditions
isDirichlet = false(Nf,1);
isNeumann = false(Nf,1);

% Neumann condition by default
isBoundary = any(G.faces.neighbors == 0,2);
isDirichlet(isBoundary) = 1;

% Find Dirichlet boundaries
if ~isempty(opt.bc)
    pressureBC = strcmpi(bc.type,'pressure');
    isDirichlet(bc.face(pressureBC)) = 1;
    isNeumann(isDirichlet) = 0;
end

%% Find subface topology for those cells that needs to be accounted for
[cno, nno, afno, fno, hfno,bfno] = createSubcellTopology(G);

[~,pos] = gridCellNodes(G,1:Nc);
% 
nCellNodes = diff(pos);
% 
% fullCell = rldecode(1:Nc,diff(pos),2)';
% 
% % Nodes here may be different from G.nodes if we ever introduce the notion
% % of calculations only on subdomains
% nodes = unique(c2n);
% nN = numel(nodes);
% n2c = cell(nN,1);
% 
% [sc2n,si] = sort(c2n); dsc2n = diff(sc2n);  d = [0; find(dsc2n > 0); numel(si)];
% for iter1 = 1 : Nn
%     n2c{iter1} = fullCell(si(d(iter1)+1 : d(iter1 + 1)));
% end
% 
% clear counter sc2n si fullCell d dsc2n c2n pos

ncno = numel(cno);

isInner = ~isBoundary(fno);
isDirichlet = isDirichlet(fno);
isNeumann = isNeumann(fno);

[~,~,subcno] = unique([cno,nno],'rows');

nsubcno = max(subcno);
nhfno = numel(unique(hfno));

nFaceNodes = diff(G.faces.nodePos);

ev = zeros(Nc,1);
for iter1 = 1 : Nc
    ev(iter1) = max(abs(eig(C{iter1})));
end

%% Neighbors of active faces

% Rows in constitutive relation:
%   contribution to x-equation from normal_x
%   contribution to y-equation from normal_x
%   ..
%   contribution to x-equation from normal_y
%   contribution to y-equation from normal_y
%  etc
% Columns: gradient of x-displacement, gradient of y-displacement, etc
%

n = G.faces.normals(fno,:);

nC = cellfun(@(x) repmat(eye(Nd),1,Nd) * diag(reshape(repmat(x,Nd,1),[],1)),num2cell(bsxfun(@rdivide,n,nFaceNodes(fno)),2),'UniformOutput',0);
tmpB = cell(ncno,1);
[tmpB{:}] = deal(C{cno});
nC = cellfun(@mtimes,nC,tmpB,'UniformOutput',0);

% Positive normal points out of the cell
normalSgn = 2 * (G.faces.neighbors(fno,1) == cno) - 1;

clear tmpB 

%% Normal stress continuity
sCont = bsxfun(@times,cat(1,nC{:}),normalSgn(reshape(repmat((1:size(normalSgn,1)),Nd,1),[],1)));

%% Discretization of normal stresses for each half face
[~,uniqueHfno] = unique(hfno);

dVal = cat(1,nC{uniqueHfno});
dRow = repmat(reshape(bsxfun(@minus,repmat(hfno(uniqueHfno)*Nd,1,Nd),fliplr(0:(Nd-1)))',[],1),1,Nd^2);
dCol = bsxfun(@minus,repmat(Nd^2*subcno(uniqueHfno),1,Nd^2),fliplr(0 : (Nd^2 - 1)) );
dCol = dCol(reshape(repmat((1:size(dCol,1)),Nd,1),[],1),:);

hook = sparse(dRow,dCol,dVal,Nd * numel(uniqueHfno),nsubcno*Nd^2);


clear dVal dRow dCol r c v

%% Continuity of displacements
if opt.weakCont
    if Nd == 2
        eta_set = [.5-sqrt(1/3)/2, .5+sqrt(1/3)/2];
        neta = 2;
    else % 3D
        eta_set = [0.5-sqrt(1/3)/2, 0.5+sqrt(1/3)/2, 0.5, 0.5;
            0, 0, -sqrt(1/3)/2, +sqrt(1/3)/2]; % Quadrature points not completely accurate in 3D (should not matter much).
        neta = 4;
    end
else
    eta_set = ones(Nd - 1,1) * opt.eta;
    neta = 1;
end

% Coordinates and differences
fc = G.faces.centroids(fno,:);
cc = G.cells.centroids(cno,:);
nc = G.nodes.coords(nno,:);
fc_cc = fc - cc;
fc_nc = fc - nc;

fn = bsxfun(@rdivide,G.faces.normals(fno,:),G.faces.areas(fno));

% Estimate for local grid size
locDX = accumarray(nno,sqrt(sum((cc - nc).^2,2)),[],@max);
locCmax = accumarray(nno,ev(cno),[],@max);

% Inner faces - one or more continuity points
% Continuity for each component of displacement for each half face
% for each continuity point (in that order)
displContPt = [];
displCC = [];

for iter1 = 1 : neta
    pt = fc_cc;
    pt(isInner,:) = pt(isInner,:) - eta_set(1,iter1) * fc_nc(isInner,:);
    ccv = zeros(numel(nno),1);
    ccv(isInner) = normalSgn(isInner);
    if size(eta_set,1) == 2
        fnXfnc = loc_cross(fn,fc - nc);
        pt(isInner,:) = pt(isInner,:) + eta_set(2,iter1) * fnXfnc(isInner,:);
    end
    if iter1 == 1
        % Use a single continuity point for Dirichlet boundaries, make this
        % the face center. That is, do nothing here.
        
        % For Dirichlet faces, the cell is always on the positive side
        ccv(isDirichlet) = 1;
    else
        pt(isDirichlet,:) = 0;
    end
    
    % Move one of the cells to the other side
    pt(isInner,:) = bsxfun(@times,normalSgn(isInner),pt(isInner,:));
    displContPt = [displContPt ; pt]; %#ok<AGROW>
    displCC = [displCC ; ccv]; %#ok<AGROW>
end

clear pt fnXfnc pt fc_nc cc nc fc_cc fn ccv ev eta_set

%% Identify and solve local systems

rInner = cell(Nn,1); cInner = rInner; vInner = rInner;
rBound = cell(Nn,1); cBound = rBound; vBound = rBound;
rDiv = cell(Nn,1); cDiv = rDiv; vDiv = rDiv;

rGrad = cell(Nn,1); cGrad = rGrad; vGrad = rGrad;
rStab = cell(Nn,1); cStab = rGrad; vStab = rStab;
rGradCorr = cell(Nn,1); cGradCorr = rGradCorr; vGradCorr = rGradCorr;

rDivBound = rDiv; cDivBound = cDiv; vDivBound = vDiv;
rGradBound = rDiv; cGradBound = cDiv; vGradBound = vDiv;
rStabBound  = cell(Nn,1); cStabBound  = rGrad; vStabBound  = rStab;
rGradCorrBound  = cell(Nn,1); cGradCorrBound = rGradCorr; vGradCorrBound  = rGradCorr;


% For convenience,
hfDir = isDirichlet(uniqueHfno);
hfNeu = isNeumann(uniqueHfno);
hf2fInd = fno(uniqueHfno);

[snno,nnoMap] = sort(nno);
dsnno = diff(snno); di = [0;find(dsnno > 0);numel(nno)];

clear snno dsnno

% Map from local (sub)-cell to global coll indices
l2gcell = 0 * subcno;

flipCell = fliplr(0:(Nd-1));
flipGrad = fliplr(0:(Nd^2-1));

trace = reshape(eye(Nd),1,Nd^2);

% Todo: 
for iter1 = 1 : Nn
    hit = nnoMap(di(iter1)+1:di(iter1+1));
    ldir = isDirichlet(hit);
    lneu = isNeumann(hit);
    
    hasStress = hit;
    hasStress(ldir) = [];
    [hfStressInd,~,riFlux] = uniqueSortedSingleOut(hfno(hasStress));
    
    hasDispl = hit;
    hasDispl(lneu) = [];
    [hfDisplInd,~,riPot] = uniqueSortedSingleOut(hfno(hasDispl));
    hp2 = hasDispl;
    
    neuRows = find(hfNeu(hfStressInd));
    isDir = hfDir(hfDisplInd);
    dirRows = find(isDir);
    nNeu = numel(neuRows);
    nDir = numel(dirRows);
    nBound = nNeu + nDir;
    
    % Extra rows for displacement continuity on internal boundaries
    if neta > 1
        intPot = hfDisplInd(~isDir); % Internal cells with potential continuity
        intPot = bsxfun(@plus,repmat(intPot,1,neta-1),nhfno * (1 : (neta - 1)));
        hfDisplInd = [hfDisplInd ; intPot(:)]; %#ok<AGROW>
        if any(~isDirichlet(hasDispl))
            notDir = hp2(~isDirichlet(hasDispl));
            hpExtra = reshape(bsxfun(@plus,repmat(notDir,1,neta-1), numel(nno) * (1:(neta - 1))),[],1);
            [~,~,hfeCP] = uniqueSortedSingleOut(hfno(notDir));
            hfeCP = reshape(bsxfun(@plus,repmat(hfeCP,1,neta-1),max(hfeCP) * (0:(neta-2))),[],1);
            
            % Recover face index for extra continuity points
            hpe2hf = mod(hpExtra,ncno);
            hpe2hf(hpe2hf == 0) = ncno;
            
            hp2 = [hp2 ; hpExtra]; %#ok<AGROW>
        else
            hpe2hf = [];
            hfeCP = [];
        end
    end
    displRows = reshape(bsxfun(@minus,repmat(Nd*hfDisplInd,1,Nd),flipCell)',[],1);

    locCells = uniqueSortedSingleOut(subcno(hit));
    globCells = uniqueSortedSingleOut(cno(hit));
    
    nlc = numel(locCells);
    nstress = Nd * numel(hfStressInd);
    ndispl = numel(displRows);
    nGrad = nlc * Nd * Nd;

    l2gcell(locCells) = 1:nlc;

    % Trivial matrix blocks
    m12 = zeros(nstress,nlc * Nd);
    m31 = zeros(nlc * Nd,nGrad);
    m32 = eye(nlc * Nd);

    % Matrix block for stress continuity
    if any(hasStress)
        srgrad = reshape(bsxfun(@minus,repmat(Nd*hasStress,1,Nd),flipCell)',[],1);
        fr = repmat(reshape(bsxfun(@minus,repmat(riFlux * Nd,1,Nd),flipCell)',[],1),1,Nd^2);
        ci = l2gcell(subcno(hasStress));
        cr = bsxfun(@minus,repmat(ci*Nd^2,1,Nd^2),flipGrad);
        fc = cr(reshape(ones(Nd,1) * (1:numel(hasStress)),[],1),:);
    else
        fr = [];
        fc = [];
        srgrad = [];
    end
    m11 = full(sparse(fr,fc,sCont(srgrad,:)/(locDX(iter1)*locCmax(iter1)),Nd * numel(hfStressInd),nlc * Nd^2));
    
    % Matrix block for weak displacement continuity
    if any(hasDispl)
        ri = [riPot; max(riPot) + reshape(hfeCP,[],1)];
        tr = bsxfun(@minus,repmat(ri * Nd,1,Nd),flipCell);
        dr = tr(:,ones(Nd,1) * (1:Nd));
        ci = l2gcell(subcno([hasDispl ; hpe2hf]));
        dc = bsxfun(@minus,repmat(ci*Nd^2,1,Nd^2),flipGrad);
        ricc = bsxfun(@minus,repmat(Nd*ri,1,Nd),flipCell);
        cicc = bsxfun(@minus,repmat(Nd*ci,1,Nd),flipCell);
    else
        dc = [];
        dr = [];
        ri = 0;
        ricc = [];
        cicc = [];
    end
    
    % Contribution from sub-cell gradients
    m21 = full(sparse(dr,dc,repmat(displContPt(hp2,:),1,Nd),Nd*max(ri),nlc * Nd^2)) / locDX(iter1);
    % Contribution from cell centers
    m22 = full(sparse(ricc,cicc,repmat(displCC(hp2),1,Nd),Nd*max(ri),nlc * Nd));
    
    % The full local system
    loc_sys = [m11 m12 ; m21 m22 ; m31 m32];
    
    % Right hand side
    nFaceCond = nstress + ndispl;
    loc_rhs = [zeros(nFaceCond,nlc * Nd); eye(nlc * Nd)];
    
    if any(hasStress)
        r = bsxfun(@minus,repmat(Nd*riFlux,1,Nd),flipCell);
        c = repmat(l2gcell(subcno(hasStress)),1,Nd);
        v = bsxfun(@times,n(hasStress,:),-(normalSgn(hasStress)./nFaceNodes(fno(hasStress))))/(locCmax(iter1));% * locDX(iter1));
        a = full(sparse(r,c,v,nstress,nlc));
        loc_rhs_stab = [full(sparse(r,c,v,nstress,nlc)); zeros(ndispl + nlc * Nd,nlc)] ;
        
    else 
        loc_rhs_stab = zeros(ndispl + nlc * Nd,nlc);
    end
    
    if iter1 == 2
%         a
    end
    
    %%% Solution of the local system
    
    if nDir > 0
        rd = reshape(bsxfun(@minus,repmat(Nd*dirRows,1,Nd),flipCell)',[],1);
        constraints =[1:nstress, nstress + rd', nstress + ndispl + (1:nlc * Nd)];
        minvars = setdiff(1:(nstress+ndispl),constraints); 
    else
        minvars = nstress + (1:ndispl);
        constraints = [1:nstress, nstress + ndispl + (1:nlc * Nd)];
    end
    
    local_mins = loc_sys(minvars,:);
    
    local_square = local_mins' * local_mins;
    rhs_square = local_mins' * loc_rhs(minvars,:);
    rhs_square_stab = local_mins' * loc_rhs_stab(minvars,:);
    
    local_constraints = loc_sys(constraints,:);
    
    loc_sys2 = [local_square , local_constraints'; ...
        local_constraints , sparse(numel(constraints),numel(constraints))];
    
    loc_ss = full(loc_sys(1:nstress,:));
    [U,S,V] = svd(loc_ss');
    
    dS = diag(S);
    threshold = 10 * eps * max(dS);
    nzS = find(dS > threshold);
    ns = numel(nzS);
    if ns == nstress
        % System is not overdetermined
        loc_rhs2 = [rhs_square ; loc_rhs(constraints,:)];
        local_solve = loc_sys2 \ loc_rhs2;
        
        loc_rhs2_stab = [rhs_square_stab ; loc_rhs_stab(constraints,:)];
        local_solve_stab = loc_sys2 \ loc_rhs2_stab;
        
        % Will miss this part when considering Neumann
%     elseif nlc < Nd % ??
%         S2 = S;
%         S2((ns+1):end,(ns+1):end) = eye(nstress-numel(nzS));
%         loc_ss2 = U * S2 * V';
%         loc_sys3 = loc_sys;
%         loc_sys3(1:nstress,:) = loc_ss2';
%         local_solve = loc_sys3 \ loc_rhs;
    else
        loc_ss2 = V(:,nzS)' * loc_ss;
        loc_con2 = [loc_ss2 ; local_constraints(nstress+1:end,:)];
        loc_sys3 = [local_square , loc_con2' ; ...
            loc_con2, sparse(size(loc_con2,1),size(loc_con2,1))];
        loc_rhs2 = [rhs_square ; V(:,nzS)' * loc_rhs(constraints(1:nstress),:); loc_rhs(constraints(nstress+1:end),:)];
        local_solve = loc_sys3 \ loc_rhs2;
        
        loc_rhs2_stab = [rhs_square_stab ; V(:,nzS)' * loc_rhs_stab(constraints(1:nstress),:); loc_rhs_stab(constraints(nstress+1:end),:)];
        local_solve_stab = loc_sys3 \ loc_rhs2_stab;
    end
    
%     if iter1 == 5
% %         local_solve_stab
% loc_rhs2_stab
%     end
    
    cellGradInd = reshape(bsxfun(@minus,repmat(Nd^2*locCells,1,Nd^2),flipGrad)',[],1);
    col = repmat(reshape(bsxfun(@minus,repmat(globCells * Nd,1,Nd),flipCell)',1,[]),nGrad,1);
    
    rInner{iter1} = reshape(repmat(cellGradInd,1,nlc * Nd),[],1);
    cInner{iter1} = reshape(col,[],1);
    vInner{iter1} = reshape(local_solve(1:nGrad,:),[],1) / locDX(iter1);

    rGrad{iter1} = reshape(repmat(cellGradInd,1,nlc),[],1);
    cGrad{iter1} = reshape(repmat(globCells',nGrad,1),[],1);
    vGrad{iter1} = reshape(local_solve_stab(1:nGrad,:),[],1) / locDX(iter1);
    
    % Note that scaling with locDX is included in div
    divOp = bsxfun(@times,kron(eye(nlc),trace),G.cells.volumes(globCells)./( nCellNodes(globCells) ))/ locDX(iter1);
    
    rDiv{iter1} = reshape(repmat(globCells,1,nlc * Nd),[],1);
    cDiv{iter1} = reshape(repmat(reshape(bsxfun(@minus,repmat(Nd*globCells,1,Nd),flipCell)',1,Nd*nlc),nlc,1),[],1);
    vDiv{iter1} = reshape(divOp * local_solve(1:nGrad,:),[],1);

    rStab{iter1} = reshape(repmat(globCells ,1,nlc),[],1);
    cStab{iter1} = reshape(repmat(globCells',nlc,1),[],1);
    vStab{iter1} = reshape(-divOp * local_solve_stab(1:nGrad,:),[],1);
    
    % Finally, the treatment of boundary conditions
    if nBound > 0
        
        fIndOfNeu = hf2fInd(hfStressInd(neuRows));
        fIndOfDir = hf2fInd(hfDisplInd(dirRows));
        
        if nNeu > 0
            rn = reshape(bsxfun(@minus,repmat(Nd*neuRows,1,Nd)',flipCell),[],1);
        else
            rn = [];
        end
        cn = 1:nNeu * Nd;
        vn = reshape(repmat((G.faces.areas(fIndOfNeu)./nFaceNodes(fIndOfNeu))',1,Nd),[],1);
        neuRHS = full(sparse(rn,cn,vn,nstress,nBound * Nd))/(locDX(iter1)*locCmax(iter1));
        
        if nDir > 0
            rd = reshape(bsxfun(@minus,repmat(Nd*dirRows,1,Nd),flipCell)',[],1);
        else
            rd = [];
        end
        cd = nNeu * Nd + (1 : nDir * Nd);
        dirRHS = full(sparse(rd,cd,1,ndispl,nBound*Nd));
        
        boundRHS = [neuRHS ; dirRHS ; zeros(nlc * Nd,nBound * Nd)];
        
        boundRHS_square = local_mins' * boundRHS(minvars,:);
        if ns == nstress
            boundRHS2 = [boundRHS_square ; boundRHS(constraints,:)];
            boundSol = (loc_sys2 \ boundRHS2) / locDX(iter1);
        else
            boundRHS2 = [boundRHS_square ; V(:,nzS)' * boundRHS(constraints(1:nstress),:); boundRHS(constraints(nstress+1:end),:)];
            boundSol = (loc_sys3 \ boundRHS2) / locDX(iter1); 
        end
        
        cb = repmat(reshape(bsxfun(@minus,repmat([fIndOfNeu ; fIndOfDir] * Nd,1,Nd),flipCell)',1,[]),nGrad,1);

        rBound{iter1} = reshape(repmat(cellGradInd,nBound*Nd,1),[],1);
        cBound{iter1} = reshape(cb,[],1);
        vBound{iter1} = reshape(boundSol(1:nGrad,:),[],1);
    end
    
    l2gcell(locCells) = 0;
end

%% Build full stress discretization
cc2subGrad = sparse(cat(1,rInner{:}),cat(1,cInner{:}),cat(1,vInner{:}));
hStress = hook * cc2subGrad;

hfi = bsxfun(@minus,repmat(hfno(uniqueHfno) * Nd,1,Nd),flipCell);
fi = bsxfun(@minus,repmat(fno(uniqueHfno) * Nd,1,Nd),flipCell);

hf2f = sparse(fi,hfi,1);

stress = hf2f * hStress;


gradPsubGrad = sparse(cat(1,rGrad{:}),cat(1,cGrad{:}),cat(1,vGrad{:}));
gradP = hf2f * hook * gradPsubGrad;


v = cat(1,nC{:});
r = repmat(reshape(bsxfun(@minus,repmat(hfno*Nd,1,Nd),fliplr(0:(Nd-1)))',[],1),1,Nd^2);
c = bsxfun(@minus,repmat(Nd^2*subcno,1,Nd^2),fliplr(0 : (Nd^2 - 1)) );
c = c(reshape(repmat((1:size(c,1)),Nd,1),[],1),:);

r = repmat((1:Nd^2 * nsubcno)',1,Nd^2);

hook2 = sparse(r,c,v);

c2 = reshape(bsxfun(@minus,repmat(cno *Nd,1,Nd),fliplr(0:(Nd-1)))',[],1);
div2 = sparse(c2,(1:ncno * Nd)',reshape(repmat(normalSgn,1,Nd)',[],1),Nc*Nd,ncno*Nd);

gradP = div2 * hook2 * gradPsubGrad;


r = bsxfun(@minus,repmat(Nd * hfno(~isDirichlet(uniqueHfno)),1,Nd),fliplr(0:(Nd-1)));
c = repmat(cno(~isDirichlet(uniqueHfno)),1,Nd);
v = bsxfun(@times,n(~isDirichlet(uniqueHfno),:),normalSgn(~isDirichlet(uniqueHfno))./nFaceNodes(fno(~isDirichlet(uniqueHfno))));

% 
% r = bsxfun(@minus,repmat(Nd * hfno(uniqueHfno(~hfDir)),1,Nd),fliplr(0:(Nd-1)));
% c = repmat(cno(uniqueHfno(~hfDir)),1,Nd);
% v = bsxfun(@times,n(uniqueHfno(~hfDir),:),normalSgn(uniqueHfno(~hfDir))./nFaceNodes(fno(uniqueHfno)));
% gpCorr = gradP + hf2f * sparse(r,c,v,size(hf2f,2),Nc);
% gp = hf2f*gpCorr;
% gradP = gradP + hf2f*sparse(r,c,v,size(hf2f,2),Nc);
gpCorr = [];

assert(norm(gradP * ones(Nc,1)) < sqrt(eps))

% Discretization of boundary conditions
row = cat(1,rBound{:});
col = cat(1,cBound{:});
val = cat(1,vBound{:});

cc2subGradBound = sparse(row,col,val,Nd^2 * nsubcno,Nf*Nd);

boundHStress = hook * cc2subGradBound;



% Bound flux gives fluxes induced by boundary conditions
boundStress = hf2f * boundHStress;

%% Divergence operating on the vector equation
tmp = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2).';
cellNo = bsxfun(@minus,repmat(Nd*tmp,1,Nd),flipCell);
cf     = bsxfun(@minus,repmat(Nd*G.cells.faces(:,1),1,Nd),flipCell);
nc     = G.cells.num;
nf     = G.faces.num;
sgn = repmat(2 * (G.faces.neighbors(G.cells.faces(:,1),1)==tmp)-1,1,Nd);
div = sparse(cf,cellNo,sgn,Nd*nf,Nd*nc)';

%%
divD = sparse(cat(1,rDiv{:}),cat(1,cDiv{:}),cat(1,vDiv{:}));
stabDelta = sparse(cat(1,rStab{:}),cat(1,cStab{:}),cat(1,vStab{:}));


%% Construct and assign output
out.stress = stress;
out.div = div;
out.A = div * stress;
out.boundStress = boundStress;
out.divD = divD;
out.stabDelta = stabDelta;
out.gradP = gradP;
% out.gp = div * gp;
out.fluidStress = gpCorr;

function c = loc_cross(a,b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2),...
    a(:,3).*b(:,1)-a(:,1).*b(:,3),...
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];
