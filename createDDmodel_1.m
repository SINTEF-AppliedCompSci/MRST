function [model_out,W,indexs] = createDDmodel_1(model,nc,Graph,W_ref)

 num_edges = numedges(Graph);
 num_nodes = numnodes(Graph);
 
  vertical_internal_faces = (nc-1)*num_edges;
horizontal_internal_faces = nc*(num_edges-1);
           internal_faces = vertical_internal_faces+horizontal_internal_faces;
           
  array_edges =  table2array(Graph.Edges);
    
  if_well     = ones(num_edges,1);
  if_inj_well = ones(num_nodes,1);
     inj_well_cell = zeros(num_nodes,1);
     
     
N = model.operators.N;
T = model.operators.T;
pv= model.operators.pv;

%% Position inticators
             k = 1;
   k_new_faces = 1;
  k_faces_zero = vertical_internal_faces;
  
%%  
           W = [];
  faces_zero = []; 
  
 for i = 1:num_edges
     if if_well(i)==1
         Well_p = array_edges(i,2);
         
         % Make producer in level k
         W = verticalWell(W, model.G, model.rock, model.G.cartDims(1), model.G.cartDims(2), k, ...      
             'Type', W_ref(Well_p).type, 'Val', W_ref(Well_p).val, ... %  SAIGUP model
            'Radius', W_ref(Well_p).r(1),'Name',W_ref(Well_p).name,'Comp_i', W_ref(Well_p).compi,'sign', W_ref(Well_p).sign);
                    %'Type', 'bhp' , 'Val', 395*barsa, ... %  Egg model 

         producer_cell = W(end).cells;
         fprintf('\n\nWe are at well producer %d with well cell % d\n', Well_p,producer_cell)

         Well_I = neighbors(Graph,Well_p);
        
         for j =1:length(Well_I)
              fprintf('\n+++ We are at conection  I%d-P%d\n',Well_I(j), Well_p);
              
                  vertical_face_zero = [];
                  
                  face_index = [];
                  face1 =  1+(k-1)*(nc-1);
                  faceN = (nc-1)+(k-1)*(nc-1); 
                  
                  cell1 = 2+(k-1)*nc;
                  cellN = (nc-1)+(k-1)*nc;
                  cell_index = [ cell1 : cellN ];
                  
                  
                  
            % Ask if we need an injector   
             if if_inj_well(Well_I(j))==1  % if Yes Make injector
                 fprintf('\n+++ Create injector and conectionr for well %d\n',Well_I(j))
                 W = verticalWell(W, model.G, model.rock, 1, 1, k, ...       
                                    'Type', W_ref(Well_I(j)).type , 'Val',  W_ref(Well_I(j)).val, ... %   
                                    'Radius', W_ref(Well_I(j)).r(1),'Name',W_ref(Well_I(j)).name,'Comp_i', W_ref(Well_I(j)).compi,'sign',W_ref(Well_I(j)).sign);
                 inj_well_cell(Well_I(j)) =  W(end).cells;
                 % Adding the face number to the conection          
                 face_index = [face_index,face1]; 
                  % Declaring Injector well as created already
                 if_inj_well(Well_I(j))=0;

             else  % if No  Conect to the corresponding injector cell
                 fprintf(2,'\n+++ Create just conection for wellpair  %s-P%s\n',W_ref(Well_I(j)).name, W_ref(Well_p).name)
                 
                 face1_new =  internal_faces + k_new_faces;
                 
                 N(face1_new,:) = [inj_well_cell(Well_I(j)) , cell1];
                 T(face1_new)   =  T(1) ;
                 
                 vertical_face_zero = [vertical_face_zero, face1];
                 
                 face_index = [face_index,face1_new]; 
                 k_new_faces = k_new_faces +1;
             end
               
             
             face_index = [face_index, face1+1 :faceN-1];
             
            % Ask if the producer was created for the first time
             if if_inj_well(Well_p)==1
                  fprintf('\n+++ Producer well number %d is created for the first time \n',Well_p)
                  face_index = [face_index,faceN];
                  if_inj_well(Well_p)=0;                 
             else
                  fprintf(2,'\n+++ Producer  well %d was already created, conection is need it \n',Well_p)
                  faceN_new =  internal_faces + k_new_faces;
                  
                  vertical_face_zero = [vertical_face_zero, faceN];
                  
                  N(faceN_new,:) = [cellN , producer_cell];
                  T(faceN_new)   =  T(1) ;
                  face_index = [face_index,faceN_new];
                  k_new_faces = k_new_faces +1;
             end
             
             if k>1
                 % Indicate the index of the previous row of horizontal faces
                 faces_zero = [faces_zero, k_faces_zero+1:k_faces_zero+nc,...'
                                           vertical_face_zero];
                 k_faces_zero =  k_faces_zero + nc;
             end

                  wellpair_cells_inx{k} = cell_index;
                  wellpair_faces_inx{k} = face_index;
                  k = k+ 1;
         end
         
        %Here we declared wich edges in the graph related the CURENT
        %PRODUCER have been  considered. With that we don't create the same
        %producer well several times.
         
         if_well( outedges(Graph,Well_p) ) = 0; 
     end
     
     
 end
 
 %% This index are for egg model
 %% TODO: improve this
indexs.faces = wellpair_faces_inx;
indexs.cells = wellpair_cells_inx;

W = W([2 3 4 5  8 6 10 12 1 7 9 11]);

%%
 
 
 model.operators =  setupOperatorsTPFA(model.G, model.rock,'neighbors',N,'trans',T,'porv',pv);
 
 % closing horizontal transmisibilities;
  model.operators.T(faces_zero) = 0;
  model_out = model;
 %% Plotting the data driven model
  figure
  plotGrid(model.G,'FaceColor', 'none')
  plotWell(model.G,W)
  view(0,0)
 
 