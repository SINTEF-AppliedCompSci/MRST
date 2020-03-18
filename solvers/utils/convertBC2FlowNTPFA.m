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
bc_nfvm.face = boundaryFaces(G);
bc_nfvm.type = repmat({'flux'}, [numel(bc_nfvm.face), 1]);
bc_nfvm.value = repmat({@(x)0}, [numel(bc_nfvm.face), 1]);

% Fill the other data: Only Neumann 1 or hom Dirichlet pressure

for i = 1:numel(bc.face)
    bcf = bc.face(i);
    j = find(bc.face(i) == bc_nfvm.face);
    bc_nfvm.type(j) = bc.type(i);
    if strcmpi(bc.type(i), 'flux')
        %[G.faces.areas(i) G.faces.areas(j)]
        %bc_nfvm.value(j) = {@(x) -1*G.faces.areas(i)}; % Neumann
        %bc_nfvm.value(j) = {@(x) bc.value(i)}; % Neumann
        bc_nfvm.value(j) = {@(x) 1 * G.faces.areas(i)}; % Neumann
    elseif strcmpi(bc.type(i), 'pressure')
        %bc_nfvm.value(j) = {@(x) bc.value(i)}; % Dirichlet pressure
        bc_nfvm.value(j) = {@(x) 1}; % Dirichlet pressure
    else
        warning('unknown bc.type')
        bc.type(i)
        keyboard
    end
end

end
