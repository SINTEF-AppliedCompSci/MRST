function G_global=fracturematrixNNC3D(G_global,tol)
% FRACTURE MATRIX NNC GENERATOR - This function takes in a global grid
% which contains G_global.Matrix and G_global.FracGrid. Fracture matrix
% intersections are detected and saved in G_global.nnc.

% METHODOLOGY - Create a loop which goes through every single fracture in
% the domain. For each fracture, submit the fracture grid and matrix grid
% to function fracmatnnc which will generate the necessary NNCs and provide
% output.

% opt=struct;
% opt=merge_options(opt,varargin{:});
% disp(['Running ',opt.type,' preprocessor...']);
t1=clock;
Gmat=G_global.Matrix;
FracGrid=G_global.FracGrid;

% Construct nnc field if not done before
if ~isfield(G_global.nnc,'cells')
    G_global.nnc=struct;
    G_global.nnc.cells=[];
    G_global.nnc.T=[];
    G_global.nnc.type={};
    G_global.nnc.area=[];
end

if isfield(G_global.rock,'poro')
    pv = poreVolume(G_global,G_global.rock);
else
    pv = G_global.cells.volumes;
end

disp('Processing Fracture-Matrix NNCs...');
for i=1:numel(fieldnames(FracGrid)) % loop through every fracture
    fieldname=['Frac',num2str(i)];
    if FracGrid.(fieldname).matrixnnc
        disp(['Records show ',fieldname,' has been processed previously for fracture-matrix NNCs.']);
        continue; % skip if fracture-matrix nnc has been processed before
    end
    disp(['Processing ',fieldname,'...']);
    F=FracGrid.(fieldname);
    
    % Detect NNCs
    [cellsi,CIi,~,typei]=fracmatnnc(F,Gmat,tol);
    
    % Compute NNC Transmissibilities
    w1 = pv(cellsi(:,1))./G_global.rock.perm(cellsi(:,1));
    w2 = pv(cellsi(:,2))./G_global.rock.perm(cellsi(:,2));
    wt = pv(cellsi(:,1))+pv(cellsi(:,2));
    % No weighting by PV:
    % w1 = 1./G.rock.perm(G.nnc.cells(:,1));
    % w2 = 1./G.rock.perm(G.nnc.cells(:,2));
    % wt = 1;
    Ti = CIi.*(wt./(w1+w2));
    
    % Adjust Ti for surface area. Ti was previously calculated using
    % frac-mat intersection area. For fractures within a cell, surface area
    % is 2*intersection. For fractures on a cell face, surface area is
    % equal to intersection area. Refer to Tene et al. (2017) pEDFM.
    Ti = adjustT(Ti, cellsi);
    
    % Append data to existing nnc information
    G_global.nnc.cells=[G_global.nnc.cells;cellsi];
    G_global.nnc.type=[G_global.nnc.type;typei];
    G_global.nnc.T=[G_global.nnc.T;Ti];  
    G_global.nnc.area=[G_global.nnc.area;F.matrix_connection.area];
    
    % Change tag to note that Frac_i has been processed
    G_global.FracGrid.(fieldname).matrixnnc=true;
end

t2=clock;
e=etime(t2,t1);
disp(['Processing Fracture-Matrix NNCs completed in ',num2str(e),' seconds!']);
    
end




