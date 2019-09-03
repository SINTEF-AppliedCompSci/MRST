function bc_nfvm = convertBC2FlowNTPFA(G, bc)

%{

Convert a std mrst bc structure to bc_nfvm for the _original_
FlowNTPFA.m. 

The bc must be such that we can loop over all faces:

for i_face=1:G.faces.num
...
ind=find(bc.face==i_face,1);
if(strcmpi(bc.type{ind},'pressure'))

That means we typically have sth like

bc_nfvm.face=boundaryFaces(G);
bc_nfvm.type=repmat({'flux'},[numel(bc_nfvm.face),1]);
bc_nfvm.value=repmat({@(x)0},[numel(bc_nfvm.face),1]);

%}

% Setup hom Neumann first. 
    bc_nfvm.face=boundaryFaces(G);
    bc_nfvm.type=repmat({'flux'},[numel(bc_nfvm.face),1]);
    bc_nfvm.value=repmat({@(x)0},[numel(bc_nfvm.face),1]);

    % Fill the other data: Only Neumann 1 or hom Dirichlet pressure 
    
    for i = 1:numel(bc.face)
        bcf = bc.face(i);
        j = find(bc.face(i) == bc_nfvm.face);
        bc_nfvm.type(j) = bc.type(i);
        if strcmpi(bc.type(i), 'flux')
            %[G.faces.areas(i) G.faces.areas(j)]
            %bc_nfvm.value(j) = {@(x) -1*G.faces.areas(i)}; % Neumann
            bc_nfvm.value(j) = {@(x) bc.value(i)}; % Neumann
        else
            bc_nfvm.value(j) = {@(x) 0}; % hom Dirichlet pressure
        end
    end

end

    
    
    
% % Convert a std mrst bc structure to bc_nfvm.
% % Fill the bc_nfvm by setting default homogeneous Neumann conditions,
% % then copy from bc_std. We should set BC explicitly on all boundary
% % faces, which means that bc_nfvm members should be of size
% % numel(boundaryFaces(G)).
% bf = boundaryFaces(G);
% nf = numel(bf);
% bc_nfvm.face = zeros(nf, 1);
% bc_nfvm.type = repmat({'flux'}, [nf, 1]);
% bc_nfvm.value = repmat({@(x) 0},[nf, 1]);

% for i = 1:nf
%     faceno = bc.face(i);
%     bc_nfvm.face(faceno) = faceno;
%     bc_nfvm.type(faceno) = bc.type(i);

%     % Can't do this with @(x) since the value is not evaluated
%     %bc_nfvm.value(faceno) = {@(x) bc_std.value(i)};
    
%     if strcmpi(bc.type(i), 'flux')
%         bc_nfvm.value(faceno) = {@(x) apa}; % flux value is 1
%     else
%         bc_nfvm.value(faceno) = {@(x) 0}; % pressure value is 0
%         disp([i,faceno,bc_nfvm.type(faceno),bc_nfvm.value(faceno)])
%     end
% end

% end
