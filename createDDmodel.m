function [model_out,W,indexs] = createDDmodel(model,nc,Graph,W_ref,levels)

 
 num_edges = numedges(Graph);
 num_nodes = numnodes(Graph);
 
  vertical_internal_faces = (nc-1)*num_edges*levels;
horizontal_internal_faces = nc*(num_edges*levels-1);
           internal_faces = vertical_internal_faces+horizontal_internal_faces;
           
  array_edges = table2array(Graph.Edges);
  
  
if_well     = ones(num_edges,1);
if_inj_well = ones(num_nodes,1);
     inj_well_cell = zeros(levels,num_nodes);
     
     
N = model.operators.N;
T = model.operators.T;
pv= model.operators.pv;

%% Position inticators
             k_levels = 1:levels; % number of cells for each well
   k_new_faces = 1;
  k_faces_zero = vertical_internal_faces; % Index where to start counting horizontal faces
 
  array_edges_2 = [];
%%  
           ref_order = [];
           W = [];
  faces_zero = []; 
 for i = 1:num_edges
     if if_well(i)==1
         Well_p = array_edges(i,2); %Taking producer on row i in the list of eadges
         
         % Make producer in level k
         % TODO ADD levels here
         W = verticalWell(W, model.G, model.rock, model.G.cartDims(1), model.G.cartDims(2), k_levels, ...      
              'Type', W_ref(Well_p).type, 'Val', W_ref(Well_p).val, ... %  SAIGUP model
            'Radius', W_ref(Well_p).r(1),'Name',W_ref(Well_p).name,'Comp_i', W_ref(Well_p).compi,'sign', W_ref(Well_p).sign);
                    %'Type', 'bhp' , 'Val', 395*barsa, ... %  Egg model 
         ref_order = [ref_order;Well_p];

         producer_cell = W(end).cells;
         fprintf('\n\nWe are at well producer %d - %s  with well cell [ ', Well_p, W_ref(Well_p).name)
             for i_prod_cell =  1:length(producer_cell)
                fprintf('%i ',producer_cell(i_prod_cell))
             end
             fprintf(']\n')

         Well_I = neighbors(Graph,Well_p);
        
         for j =1:length(Well_I)
              fprintf('\n+++ We are at conection  %s - %s \n',W_ref(Well_I(j)).name , W_ref(Well_p).name);
              
                  vertical_face_zero = [];
                  
                  face_index = zeros(levels,nc-1); % Would be the vector containing the index of a articular conection
                  face1 =  1+(k_levels'-1)*(nc-1); 
                  faceN = (nc-1)+(k_levels'-1)*(nc-1); 
                  
                  cell1 = 2+(k_levels'-1)*nc;
                  cellN = (nc-1)+(k_levels'-1)*nc;
                  
                  for ii_conx =  1 : levels
                    cell_index(ii_conx,:) = [ cell1(ii_conx) : cellN(ii_conx) ];
                  end
                  
                  
                  
            % Ask if we need an injector   
             if if_inj_well(Well_I(j))==1  % if Yes Make injector
                 fprintf('\n+++ Create injector and conectionr for well %s  - %s\n', W_ref(Well_I(j)).name,W_ref(Well_p).name)
                 fprintf('\n+++ At level %d - %d \n', k_levels(1),k_levels(2))

                 % TODO ADD levels here
                 W = verticalWell(W, model.G, model.rock, 1, 1, k_levels, ...       
                                   'Type', W_ref(Well_I(j)).type , 'Val',  W_ref(Well_I(j)).val, ... %   
                                    'Radius', W_ref(Well_I(j)).r(1),'Name',W_ref(Well_I(j)).name,'Comp_i', W_ref(Well_I(j)).compi,'sign',W_ref(Well_I(j)).sign);
                 ref_order = [ref_order;  Well_I(j)];

                 inj_well_cell(:,Well_I(j)) =  W(end).cells;
                 % Adding the face number to the conection          
                 face_index(:,1) = face1; 
                  % Declaring Injector well as created already
                 if_inj_well(Well_I(j))=0;

             else  % if No  Conect to the corresponding injector cell
                 fprintf(2,'\n+++ Create just conection for wellpair  %s - %s\n', W_ref(Well_I(j)).name,  W_ref(Well_p).name)
                 
                 fprintf('\n+++ At level %d - %d \n', k_levels(1),k_levels(2))

                 face1_new =  internal_faces + k_new_faces;
                 
                 for ii_conx =   1 : levels
                     N(face1_new+ii_conx-1,:) = [inj_well_cell(ii_conx,Well_I(j)) , cell1(ii_conx)];
                     T(face1_new+ii_conx-1)   =  T(1) ;
                     face_index(ii_conx,1)    = face1_new+ii_conx-1;
                 end
                 vertical_face_zero = [vertical_face_zero, face1'];
                 
                 k_new_faces = k_new_faces +levels;
             end
               
             for ii_conx=1:levels
                 face_index(ii_conx,[2:nc-2]) = [face1(ii_conx)+1 :faceN(ii_conx)-1];
             end
             
            % Ask if the producer was created for the first time
             if if_inj_well(Well_p)==1
                  fprintf('\n+++ Producer well number %s is created for the first time \n', W_ref(Well_p).name)
                  face_index(:,end) = faceN;
                  if_inj_well(Well_p)=0;                 
             else
                  fprintf(2,'\n+++ Producer  well %s was already created, conection is need it \n',  W_ref(Well_p).name)
                  faceN_new =  internal_faces + k_new_faces;
                  
                  vertical_face_zero = [vertical_face_zero, faceN'];
                  for ii_conx =   1 : levels
                    N(faceN_new+ii_conx-1 ,:) = [cellN(ii_conx) , producer_cell(ii_conx)];
                    T(faceN_new+ii_conx-1)   =  T(1) ;                   
                    face_index(ii_conx,end) = faceN_new+ii_conx-1 ;
                  end 
                  
                  k_new_faces = k_new_faces + levels;
             end
             
             if k_levels>1
                 % Indicate the index of the previous row of horizontal faces
                 
                            % previous faces indx   , horizontal face index
                            % , vertifcal faces index
                 faces_zero = [faces_zero, k_faces_zero+1:k_faces_zero+levels*nc,...' 
                                           vertical_face_zero];
                                       
                 k_faces_zero =  k_faces_zero + nc*levels;
             end
             
             for ii_conx =  1 : levels
                 wellpair_cells_inx{k_levels(ii_conx)} = cell_index(ii_conx,:);
                 wellpair_faces_inx{k_levels(ii_conx)} = face_index(ii_conx,:);
             end
                  k_levels = k_levels+ levels;
                  array_edges_2 = [array_edges_2;Well_p, Well_I(j) ];
         end
         
        %Here we declared wich edges in the graph related the CURENT
        %PRODUCER have been  considered. With that we don't create the same
        %producer well several times.
         
         if_well(outedges(Graph,Well_p) ) = 0; 
     end
     
     
 end
 
 if levels >1
                  % Indicate the index of the previous group of row of horizontal faces

                % previous faces indx   , horizontal face index
                % , vertifcal faces index
     faces_zero = [faces_zero, k_faces_zero+1:k_faces_zero+(levels-1)*nc,...' 
                               vertical_face_zero];
                           
 end
 
 
% indexs.faces = {wellpair_faces_inx{[1:30,35 36 31:34]}};
% indexs.cells = {wellpair_cells_inx{[1:30,35 36 31:34]}};
%  
%  indexs.faces = {wellpair_faces_inx{[1:36,49:54,37:48]}};
%  indexs.cells = {wellpair_cells_inx{[1:36,49:54,37:48]}};
 
%   indexs.faces = {wellpair_faces_inx{[ 1:16 ,123:124,   47:52, 17:18,...
%                                       19:28 ,125:130, 131:134, 109:112, 99:102,...
%                                       53:60 ,61:68  , 113:118, 29:44,...
%                                     103:108 ,119:122,   69:72, 73:80, 81:88,...
%                                       89:98 ,45:46 ]}};
%   indexs.cells = {wellpair_cells_inx{[ 1:16 ,123:124,   47:52, 17:18,...
%                                       19:28 ,125:130, 131:134, 109:112, 99:102,...
%                                       53:60 ,61:68  , 113:118, 29:44,...
%                                     103:108 ,119:122,   69:72, 73:80, 81:88,...
%                                       89:98 ,45:46 ]}};
 
%% Re-order the parameter index following the order in DD.wps
   [~,I_Graph_index]=sort(array_edges_2(:,1));   
  for i =1: numel(I_Graph_index)
      indx = I_Graph_index(i); 
      for k =  1 : levels
          out_indx(levels*(i-1) + k) = levels*(I_Graph_index(i)-1)+k;
      end
  end
  indexs.faces = {wellpair_faces_inx{out_indx}};
  indexs.cells = {wellpair_cells_inx{out_indx}};
 
 model.operators =  setupOperatorsTPFA(model.G, model.rock,'neighbors',N,'trans',T,'porv',pv);
 
 % closing horizontal transmisibilities;
  model.operators.T(faces_zero) = 0;
  model_out = model;
  
  [~,I]=sort(ref_order,1);
  W = W(I); %Reordering wells in their original sequence in W_ref

  
 %% Plotting the data driven model
  figure
  plotGrid(model.G,'FaceColor', 'none')
  plotWell(model.G,W)
  view(0,0)
  
 
 