classdef IMSIMmodel 
    %IMSIMMODEL Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Graph
        W
        model
    end
    
    methods
        function obj = IMSIMmodel(model,nc,Graph,W_in)
            %IMSIMMODEL Construct an instance of this class
            %   Detailed explanation goes here
            obj.model = model;
            N = model.operators.N;
            T = model.operators.T;
            pv= model.operators.pv;
            
            number_of_vertical_faces = (nc-1)*size(Graph.Edges.EndNodes,1); %Internal faces
            T(number_of_vertical_faces+1 :end) = 0;
            Graph.Nodes.Face_Indx =0*(1:12)';%TODO: check this


            nf =  nc -1;
            for i_ed =  1:size(Graph.Edges.EndNodes,1)
                
                %% Savinf cell indexs
                fprintf('Cell indexs Edge %i, %i-%i \n',i_ed,Graph.Edges.EndNodes(i_ed,1),Graph.Edges.EndNodes(i_ed,2));
                
     
                [Cell_1 , Cell_nc] = deal( 1+(i_ed-1)*nc  ,...
                                          nc+(i_ed-1)*nc) ;
                                      
                Graph.Edges.Cell_Indices{i_ed} = Cell_1+1 : Cell_nc-1;
                
                
                %% Saving internal faces indexs
                fprintf('in-Faces indexs Edge %i, %i-%i \n',i_ed,Graph.Edges.EndNodes(i_ed,1),Graph.Edges.EndNodes(i_ed,2));
                
                [Node_1,Node_nc] = deal(Graph.Edges.EndNodes(i_ed,1),...
                                       Graph.Edges.EndNodes(i_ed,2));
                % Check if the firs node was already defined and has its
                % face number index stored
                if Graph.Nodes.Face_Indx(Node_1)== 0                    
                    
                    %fprintf(2,'Cell number %i, will have Well %s \n',Cell_1 ,W_in(Node_1).name); 
                        
                    First_Index =  1+(i_ed-1)*nf;  
                    Graph.Nodes.Face_Indx(Node_1) = First_Index;
                    Graph.Nodes.Well_cells(Node_1)= Cell_1;
                else
                    
                    %Create New conectivity
                    N = [N;...
                         Graph.Nodes.Well_cells(Node_1), Cell_1+1]; % New artificial face
                    T = [T;T(1)];
                    
                    First_Index = length(T); % This the indexa a new face  
                                             % between cells
                                             % Graph.Nodes.Well_cells(Node_1)
                                             % and  Cell_1+1
                end
                
                % Internal faces number
                Internal    =  (2+(i_ed-1)*nf:(nc-2)+(i_ed-1)*nf) ;
                
                % Check if the last node was already defined and has its
                % face number index stored
                if Graph.Nodes.Face_Indx(Node_nc) ==0
                    %fprintf(2,'Cell number %i, will have Well %s \n',Cell_nc ,W_in(Node_nc).name);                     
                    Last_Index  =  (nc-1)+(i_ed-1)*nf; 
                    Graph.Nodes.Face_Indx(Node_nc) = Last_Index;
                    Graph.Nodes.Well_cells(Node_nc)= Cell_nc;
                else
                     %Create New conectivity
                    N = [N;...
                         Cell_nc-1, Graph.Nodes.Well_cells(Node_nc)]; % New artificial face
                    T = [T;T(1)];
                    
                    Last_Index = length(T); % This the index of a new face  
                                             % between cells
                                             % Graph.Nodes.Well_cells(Node_1)
                                             % and  Cell_1+1
                end
                fprintf('\n');
                Graph.Edges.Face_Indices{i_ed} = [First_Index,Internal,Last_Index];                
            end
            
            
            W = [];
            for iw = 1 : numel(W_in)
%                 W_in(iw).cells = Graph.Nodes.Well_cells(iw);
%                 
%                 %Reducing the size of well cells to 1
%                 W_in(iw).r =W_in(iw).r(1) ;
%                 W_in(iw).rR =W_in(iw).rR(1) ;
%                 W_in(iw).WI =W_in(iw).WI(1) ;
%                 W_in(iw).cstatus =W_in(iw).cstatus(1) ;
%                 W_in(iw).cell_origin =W_in(iw).cell_origin(1) ;
%                 
                
                W = verticalWell(W, model.G, model.rock, model.G.cartDims(1), model.G.cartDims(2), 1, ...      
             'Type', W_in(iw).type, 'Val', W_in(iw).val, ... %  SAIGUP model
            'Radius',W_in(iw).r(1),'Name',W_in(iw).name,'Comp_i', W_in(iw).compi,'sign', W_in(iw).sign);
        
                W(iw).cells = Graph.Nodes.Well_cells(iw);

            end
            
            
            obj.Graph =  Graph;
            obj.W     =  W;
            %obj.model.rock.perm = [0*model.rock.perm, model.rock.perm, 0*model.rock.perm];
            
            %obj.model.operators = setupOperatorsTPFA(obj.model.G, ...
            %                                         obj.model.rock);
                                                 
           % obj.model = model.removeStateFunctionGroupings();
            
            obj.model.operators =  setupOperatorsTPFA( obj.model.G,  obj.model.rock,...
                                                      'neighbors',N,...
                                                      'trans',T,...
                                                      'porv',pv);
                                                  
             %% Plotting the data driven model
              figure
              plotGrid(obj.model.G,'FaceColor', 'none')
              plotWell(obj.model.G,obj.W)
              view(0,0)
      
        end        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

