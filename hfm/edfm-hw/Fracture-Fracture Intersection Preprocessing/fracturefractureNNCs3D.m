function [G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol,varargin)
% Use this function after G=assembleGlobalGrid(G) and 
% G=computeEffectiveTrans(G) have been called out.
%
% This function takes in the global grid structure G which at this point
% should contain G.FracGrid.FracN where N represents each fracture plane.
% Additionally, fracplanes is also taken as an input. The last argument is
% the tolerance for double point floating numbers.
% The function does the following:
% 1. Detect fracture plane intersections
% 2. Establish NNCs and calculate transmissibilities for intersecting 
%    fracture grids
% 3. Save NNCs in G.nnc
% 4. Save intersection relationships in fracplanes.intersects
% 
% Returns G with additional NNCs in G.nnc
% Returns fracplanes with information about intersections with other
% fracture planes.
%
% Future improvements: Add in G.nnc.type

opt=struct('Verbose',true);
opt=merge_options(opt,varargin{:});
Verbose=opt.Verbose;

t1=clock;

if ~isfield(G.nnc,'cells')
    G.nnc=struct;
    G.nnc.cells=[];
    G.nnc.T=[];
    G.nnc.type={};
    G.nnc.area=[];
end

disp('Processing Fracture-Fracture NNCs...');
Gf=G.FracGrid;
numfrac=numel(fieldnames(Gf)); % number of fractures

% METHODOLOGY:
% Starting with Frac1, we look at Frac2,3,4...N and determine the NNCs
% between Frac1 and the other Fracs.
% For Frac2, since Frac1&2 have previously been studied, Frac2 only needs
% to be studied with Frac3,4,5...N.
% This pattern continues, in that Fraci only needs to be studied with
% Frac(i+1)...N.
% Finally, Frac(N-1) is studied with FracN and the process stops.
% So, the outermost loop goes from i=1:(N-1)
% The first inner loop goes from j=(i+1):N
%
% For every Frac(i) and Frac(j) pair, the first quick check is to see if
% they are parallel. If they are, then they do not intersect (because when
% setting up the fracture planes, it does not make sense to overlap them
% anyway). If they are not parallel, then call function fgridNNCs3D which
% will append new NNCs to G and saves intersection information in 
% fracplanes.

% loop through each fracture plane
for i=1:(numfrac-1)
    % loop through every other fracture plane
    fieldnamei=['Frac',num2str(i)];
    
    % Frac_i corresponds to fracplanes_k
    k=Gf.(fieldnamei).fracplanenumber;
    
    mcells_i=Gf.(fieldnamei).matrix_connection.cells;
    cloudlist=cell(size(mcells_i));
    
    
    for j=(i+1):numfrac
        fieldnamej=['Frac',num2str(j)];
        
        % check if fraci and fracj NNCs have been processed before.
        if ismember(j,G.FracGrid.(fieldnamei).fracgridnnc)
            disp([fieldnamei,' has been checked against ',fieldnamej,' previously.']);
            continue; % skip to next fracj
        end
        
        dispif(Verbose,['Checking ',fieldnamei,' vs ',fieldnamej,'.\n']);
        G.FracGrid.(fieldnamei).fracgridnnc(end+1)=j; % record that i and j have been checked
        G.FracGrid.(fieldnamej).fracgridnnc(end+1)=i;
                
        % Frac_j corresponds to fracplanes_l
        l=Gf.(fieldnamej).fracplanenumber;
        
        if isfield(fracplanes,'SetID')
            setID_i=fracplanes(k).SetID;
            setID_j=fracplanes(l).SetID;
        else
            setID_i=1;
            setID_j=2;
        end
        
        % If setID's are the same, then they do not intersect
        if setID_i==setID_j
            continue;
        end
        
        % check if Frac_i and Frac_j are parallel. If parallel, continue to
        % next iteration. This is a quick check to save processing time. If
        % the fracplanes are from the same fractureID, then they may be
        % parallel and the intersection check is not skipped
        if isfield(fracplanes,'fractureID')
            fractureID_i=fracplanes(k).fractureID;
            fractureID_j=fracplanes(l).fractureID;
        else
            fractureID_i=1;
            fractureID_j=2;
        end
        
        if fractureID_i~=fractureID_j && ...
                norm(cross(fracplanes(k).normal,fracplanes(l).normal))<tol
            continue;
        end
        
        % determine if two fracture planes intersect (given that they are
        % not parallel or from the same fractureID])
        [xsect,fullline,~,~]=PEBIPEBIintersection3D(fracplanes(k).points,fracplanes(l).points,tol);
        if ~xsect
            continue;
        end
        
        if isempty(cloudlist{1})
            for jj=1:length(mcells_i)
                mcell=mcells_i(jj);
                cloudlist{jj}=[mcell;findneighbours(G.Matrix,mcell)];
            end
        end
        
        % generate NNCs between Frac_i and Frac_j
        [G,fracplanes]=fgridNNCs3D(G,fracplanes,cloudlist,fullline,i,j,k,l,tol);
        
    end
end

% Save NNC type under G.nnc.type
numfracfracnncs=size(G.nnc.cells,1)-size(G.nnc.type,1);
type=cell(numfracfracnncs,1);
for i=1:numfracfracnncs
    type{i}='fracfrac';
end
G.nnc.type=[G.nnc.type;type];
G.nnc.area=[G.nnc.area;zeros(size(type))];

t2=clock;
e=etime(t2,t1);
disp(['Processing Fracture-Fracture NNCs completed in ',num2str(e),' seconds!']);

end
