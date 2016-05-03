clc; clear all; close all;
addpath('../');

add = .5;
lineWidth = 1.5;
faceAlpha = .2;
mrkSzBig = 20;
mrkSzSml = 10;
cut = 4;
addTxt = .25;
azel = [-30,-15];

fN = 6;

switch fN
    case 1
        xMax = 1.2; yMax = 1;
        G = cartGrid([1,1], [xMax, yMax]);
        G = computeGeometry(G);

        plotGrid(G, 'facealpha', faceAlpha)
        set(gcf, 'defaulttextinterpreter','latex');
        
        set(gca,'XTick',[0 xMax] );
        set(gca,'XTickLabel',{'', ''} )
        set(gca,'YTick',[0 yMax] );
        set(gca,'YTickLabel',{'', ''} )
        axis([-xMax*.2 xMax*1.2 -yMax*.2 yMax*1.2]) 
        my_xticklabels([0,xMax], {'$-h_x$' '$h_x$'}, .0, -.07)
        my_yticklabels([0,yMax], {'$-h_y$' '$h_y$'}, .04, -.1)
        axis equal
        
        xlabel('$x$'); ylabel('$y$');
        
        ps = get(gcf, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(gcf, 'paperunits', 'centimeters');
        set(gcf, 'papersize', [paperWidth paperHeight]);
        set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

        print(gcf, '-dpdf', '../../tex/thesis/fig/Example1.pdf');
        
    case 2
        xMax = 1.2; yMax = 1; zMax = 1.5;
        G = cartGrid([1,1,1], [xMax, yMax, zMax]);
        G = computeGeometry(G);

        plotGrid(G, 'facealpha', faceAlpha)
        set(gcf, 'defaulttextinterpreter','latex');
        
        set(gca,'XTick',[0 xMax] );
        set(gca,'XTickLabel',{'-h_x', 'h_x'} )
        set(gca,'YTick',[0 yMax] );
        set(gca,'YTickLabel',{'-h_y', 'h_y'} )
        set(gca,'ZTick',[0 zMax] );
        set(gca,'ZTickLabel',{'h_z', '-h_z'} )
        set(gca, 'FontName', 'Times-Roman', 'FontAngle', 'Oblique')
        axis equal
        axis([-xMax*.2 xMax*1.2 -yMax*.2 yMax*1.2 -zMax*.2 zMax*1.2]) 
        
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        
        view(azel)
        
        ps = get(gcf, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(gcf, 'paperunits', 'centimeters');
        set(gcf, 'papersize', [paperWidth paperHeight]);
        set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

        print(gcf, '-dpdf', '../../tex/thesis/fig/Example2.pdf');
        
    case 3
        
          
        P = [-1,-1; 1,-1; -1,1];
        G = triangleGrid(P);
        G = computeGeometry(G);
        plotGrid(G, 'facealpha', faceAlpha)
        set(gcf, 'defaulttextinterpreter','latex');
        
        set(gca,'XTick',[-1 1] );
        set(gca,'YTick',[-1 1] );
        axis equal
        axis([-1.5, 1.5, -1.5, 1.5]) 
        
        xlabel('$x$'); ylabel('$y$');

        V = G.nodes.coords;
        hold on;
        plot(V(:,1), V(:,2), 'ok', 'MarkerFaceColor', 'r')
        text(-1.2, -1.2, '$V_1$')
        text(1.1,  -1.2, '$V_2$')
        text(-1.2,  1.2, '$V_3$')
        text(-1.3, 0, '$E_3$')
        text(0, -1.2, '$E_1$')
        text(.1, .1, '$E_2$')
        
                ps = get(gcf, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(gcf, 'paperunits', 'centimeters');
        set(gcf, 'papersize', [paperWidth paperHeight]);
        set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

        print(gcf, '-dpdf', '../../tex/thesis/fig/Kt.pdf');
        
    case 4
        
        P = [-1, -1, -1; 1 -1 -1; 1 1 -1; -1 1 -1;
             -1, -1, 1 ; 1 -1  1; 1 1  1; -1 1  1];
         
        T = delaunay(P);
        G = tetrahedralGrid(P,T);
        G.nodes.coords(:,3) = -G.nodes.coords(:,3);
        G = mrstGridWithFullMappings(G);
        
        plotGrid(G,1, 'facealpha', faceAlpha)
        set(gcf, 'defaulttextinterpreter','latex');
        
        set(gca,'XTick',[-1 1] );
        set(gca,'YTick',[-1 1] );
        set(gca,'ZTick',[-1 1] );
        set(gca, 'ZTickLabels', {1, -1});
        axis equal
        axis([-1.5, 1.5, -1.5, 1.5 -1.5 1.5]) 
        
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');
        
        hold on
        nodeNum = G.cells.nodePos(1):G.cells.nodePos(2)-1;
        nodes = G.cells.nodes(nodeNum);
        V = G.nodes.coords(nodes,:);
        plot3(V(:,1), V(:,2), V(:,3), 'ok', 'MarkerFaceColor', 'r')
        text(-1.2, -1.2, 1.3,  '$V_1$')
        text(1.2, -1.2, 1.2, '$V_2$')
        text(-1.2, 1.2, 1.2, '$V_3$')
        text(-1.2, -1.2, -1.2, '$V_4$')
        text(-.2, -1.3, .9, '$E_1$')
        text(.2, .2, 1.2,'$E_2$')
        text(-1.3, -.2, .9, '$E_3$')
        text(-1, -1.2, .1, '$E_4$')
        text(.1, -1.2, -.2, '$E_5$')
        text(-1.2, .2, -.2, '$E_6$')
        
        view(azel)
        
        ps = get(gcf, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(gcf, 'paperunits', 'centimeters');
        set(gcf, 'papersize', [paperWidth paperHeight]);
        set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

        print(gcf, '-dpdf', '../../tex/thesis/fig/KT.pdf');
        
        
    case 5
        
        xMax = 2; yMax = 2;
        G = cartGrid([1,1], [xMax, yMax]);
        X = G.nodes.coords;
        X = bsxfun(@minus, X,[1,1]);
        G.nodes.coords = X;
        G = computeGeometry(G);
        plotGrid(G, 'facealpha', faceAlpha)
        set(gcf, 'defaulttextinterpreter', 'latex')
        
        set(gca,'XTick',[-1 1] );
        set(gca,'YTick',[-1 1] );
        axis equal
        axis([-1.5, 1.5, -1.5, 1.5]) 
        
        xlabel('$x$'); ylabel('$y$');

        V = G.nodes.coords;
        hold on;
        plot(V(:,1), V(:,2), 'ok', 'MarkerFaceColor', 'r')
        text(-1.2, -1.2, '$V_1$')
        text(1.1,  -1.2, '$V_2$')
        text(1.1, 1.2, '$V_3$')
        text(-1.2,  1.2, '$V_4$')
        text(-1.3, 0, '$E_1$')
        text(1.1, 0, '$E_2$')
        text(0, -1.2, '$E_3$')
        text(0, 1.2, '$E_4$')
        
                ps = get(gcf, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(gcf, 'paperunits', 'centimeters');
        set(gcf, 'papersize', [paperWidth paperHeight]);
        set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

        print(gcf, '-dpdf', '../../tex/thesis/fig/Ks.pdf');
        
            
    case 6
        
        xMax = 2; yMax = 2; zMax = 2;
        G = cartGrid([1,1,1], [xMax, yMax, zMax]);
        X = G.nodes.coords;
        X = bsxfun(@minus, X,[1,1,1]);
        G.nodes.coords = X;
        G = computeGeometry(G);
        plotGrid(G, 'facealpha', faceAlpha)
        set(gcf, 'defaulttextinterpreter','latex');
        
        set(gca,'XTick',[-1 1] );
        set(gca,'YTick',[-1 1] );
        set(gca,'ZTick',[-1 1] );
        axis equal
        axis([-1.5, 1.5, -1.5, 1.5 -1.5 1.5]) 
        
        xlabel('$x$'); ylabel('$y$'); zlabel('$z$');

        V = G.nodes.coords;
        hold on;
        plot3(V(:,1), V(:,2), V(:,3), 'ok', 'MarkerFaceColor', 'r')
        text(-1.2,-1.2, 1.2, '$V_1$')
        text( 1.1,-1.2, 1.2, '$V_2$')
        text(  .8,  .8, 1.2, '$V_3$')
        text(-1.3, 1.3, 1.2, '$V_4$')
        text(- .8,- .8,-1.2, '$V_5$')
        text( 1.1,-1.2,-1.2, '$V_6$')
        text( 1.2, 1.2,-1.2, '$V_7$')
        text(-1.3, 1.3,-1.2, '$V_8$')
        text(   0,-1.1, 1.2, '$E_1$')
        text(  .9,  .0, 1.2, '$E_2$')
        text(   0, 1.0, 1.2, '$E_3$')
        text(-1.2,   0, 1.2, '$E_4$')
        text(- .9,- .9,  .0, '$E_5$')
        text( 1.1,- .9,  .0, '$E_6$')
        text( 1.1, 1.0,  .0, '$E_7$')
        text(- .8, 1.1,  .0, '$E_8$')
        text(   0,- .9,-1.2, '$E_{9}$')
        text( 1.1,  .0,-1.2, '$E_{10}$')
        text(   0, 1.1,-1.2, '$E_{11}$')
        text(-1.0,   0,-1.2, '$E_{12}$')
        
        ps = get(gcf, 'Position');
        ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
        paperWidth = 10;
        paperHeight = paperWidth*ratio - cut;
        set(gcf, 'paperunits', 'centimeters');
        set(gcf, 'papersize', [paperWidth paperHeight]);
        set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);

        print(gcf, '-dpdf', '../../tex/thesis/fig/KC.pdf');
        
end