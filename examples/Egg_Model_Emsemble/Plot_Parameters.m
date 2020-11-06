

n_edges = numedges(Graph);
for i = 1 :n_edges
    edges_name{i} = num2str(Graph.Edges.EndNodes(i,:));
end
val_opt = control2value(p_opt,parameters);

figure('Name','Parameters values','NumberTitle','off')

  
subplot(1,3,1)

barh(categorical(edges_name) ,[val{1};val_opt{1}'])
set(gca,'YDir','reverse')
title('transmisibilities')

subplot(1,3,2)

barh(categorical(edges_name) ,[val{1};val_opt{1}'])
set(gca,'YDir','reverse')
title('Pore Volume')

subplot(1,3,3)

barh(categorical(Graph.Nodes.Well_name) ,[val{3}';val_opt{3}'])
set(gca,'YDir','reverse')
title('WI')


figure('Name','Parameters control','NumberTitle','off')

  
subplot(1,3,1)

barh(categorical(edges_name) ,[p0_fd(1:n_edges)';p_opt(1:n_edges)'])
set(gca,'YDir','reverse')
title('transmisibilities')

subplot(1,3,2)

barh(categorical(edges_name) ,[p0_fd(n_edges+1:2*n_edges)';p_opt(n_edges+1:n_edges*2)'])
set(gca,'YDir','reverse')
title('Pore Volume')

subplot(1,3,3)

barh(categorical(Graph.Nodes.Well_name) ,[p0_fd(2*n_edges+1:end)';p_opt(2*n_edges+1:end)'])
set(gca,'YDir','reverse')
title('Pore Volume')


figure
                pv = val{2};
                TT=  val{1};
                
                
                val_opt = control2value(p_opt, parameters);
                pv_opt = val_opt{2};
                TT_opt=  val_opt{1};

            n_wells  = numel(W);
            for i = 1:n_wells
                cell_number = W_ref(i).cells(1);
                XData(i) = model_ref.G.cells.centroids(cell_number,1);
                YData(i) = model_ref.G.cells.centroids(cell_number,2);
                ZData(i) = model_ref.G.cells.centroids(cell_number,3);
            end
subplot(2,2,1);

                h2= plotGrid(model_ref.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);                
                hold on, pg =  plot(Graph, 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*TT/max(TT))
                        labelnode(pg,[1:n_wells],{W.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                title("Initial Transmisibility")


subplot(2,2,2);
                h2= plotGrid(model_ref.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);                
                hold on, pg=plot(Graph,'r','XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*pv/max(pv))
                         labelnode(pg,[1:n_wells],{schedule.control(1).W.name})
                hold off;
                pg.NodeFontSize= 15;
                axis off ; 
                title("Initial Pore volume")

subplot(2,2,3);
                
                h2= plotGrid(model_ref.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);                
                hold on, pg =  plot(Graph, 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*TT_opt/max(TT))
                        labelnode(pg,[1:n_wells],{W.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                title("Optimize Transmisibility")

subplot(2,2,4);
                h2= plotGrid(model_ref.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);                
                hold on, pg=plot(Graph,'r','XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*pv_opt/max(pv))
                         labelnode(pg,[1:n_wells],{schedule.control(1).W.name})
                hold off;
                pg.NodeFontSize= 15;
                axis off ; 
                title("Optimized Pore volume")