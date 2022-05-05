function Subplots(f,figSubtitle)
    figure(f);
%     t = tiledlayout(length(figSubtitle)/2,length(figSubtitle)/2);
%     t = tiledlayout(3,2);
%     t(1) = nexttile([2 1]); title(figSubtitle{1});
%     t(2) = nexttile; title(figSubtitle{2});
%     t(3) = nexttile; title(figSubtitle{3});
%     t(4) = nexttile; % this is slider tile
%     t(5) = nexttile; title(figSubtitle{4});

%     t = tiledlayout(3,2);
%     t(1) = nexttile([2 1]); title(figSubtitle{1});
%     for i = 2 : length(figSubtitle)
%         t(i) = nexttile; title(figSubtitle{i});
%     end    

    t = tiledlayout(3,2);
    t(1) = nexttile; title(figSubtitle{1});
    t(2) = nexttile; title(figSubtitle{2});
    t(3) = nexttile; title(figSubtitle{3});
    t(4) = nexttile; title(figSubtitle{4});    
    t(5) = nexttile; title(figSubtitle{5});    
    t(6) = nexttile; title(figSubtitle{6});
   
    % if consider an inset plot in another plot
%     parentPos = t(3).Position;
%     width = parentPos(3) / 2.5; 
%     height = parentPos(4) / 3;
%     offsetX = 0.1; offsetY = 0.4;
%     left = parentPos(1) + parentPos(3) - (1 + offsetX) * width;
%     bottom = parentPos(2) + (offsetY) * height;
%     insetAx = axes('position',[left bottom width height]);
%     insetAx.Title.String = figSubtitle{end};
end