%load paper_EGG_Model_40.mat

%close all
for i = 1:numel(schedule.step.val)
     for p = 9:12
        qw_ref(i,p-8) = wellSols_ref{i}(p).qWs;
        qo_ref(i,p-8) = wellSols_ref{i}(p).qOs;
        qw_0(i,p-8) = wellSols_0{i}(p).qWs;
        qo_0(i,p-8) = wellSols_0{i}(p).qOs;
        qw_pred(i,p-8) = wellSols_opt{i}(p).qWs;
        qo_pred(i,p-8) = wellSols_opt{i}(p).qOs;
     end
     qw_T_ref(i) = sum( abs(qw_ref(i,:)));
     qo_T_ref(i) = sum( qo_ref(i,:));
     qw_T_0(i) = sum( abs(qw_0(i,:)));
     qo_T_0(i) = sum( qo_0(i,:));
     qw_T_pred(i) = sum( abs(qw_pred(i,:)));
     qo_T_pred(i) = sum( qo_pred(i,:));
 end
 tt = cumsum(schedule.step.val)/day;
 f_units =  day/stb;


%% figure paper
    
lineSize = 2;
markersize = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555555
 figure

 hs =subplot(1,2,1)
 h3 = plot(hs,tt,qw_T_ref*f_units,'.k',tt,qw_T_0*f_units,tt,qw_T_pred*f_units);
 set(h3,{'LineWidth'},{lineSize;lineSize;lineSize})

 xlabel('days'); ylabel('STB');
 %ylim([0 4000])
 title('Water production')
 T = cumsum(dt(1:48))/day; %your point goes here 
 hl = line([T(end) T(end)],get(hs,'YLim'),'Color', 150*[1 1 1]/255,'LineStyle','--','HandleVisibility','off');
 set(hl,'LineWidth',lineSize)
 %text(T(end),'Training data')
 legend('Reference','Initial model','Calibrated model','Location','SouthEast') 


 
  hs = subplot(1,2,2)
  h4 = plot( hs,tt,-qo_T_ref*f_units,'.k',tt,-qo_T_0*f_units,tt,-qo_T_pred*f_units)
  set(h4,{'LineWidth'},{lineSize;lineSize;lineSize})
  xlabel('days'); ylabel('STB'); %ylim([0 4000])
  title('Oil production')
  legend('Reference','Initial model','Calibrated model','Location','NorthEast')     
 hl = line([T(end) T(end)],get(hs,'YLim'),'Color', 150*[1 1 1]/255,'LineStyle','--','HandleVisibility','off');
 set(hl,'LineWidth',lineSize)

    
  figure  
  
  for i = 1:4
      hs  = subplot(2,2,i)
           h5 = plot(hs,tt,-qw_ref(:,i)*f_units,'.k',tt,-qw_0(:,i)*f_units,tt,-qw_pred(:,i)*f_units)
             set(h5,{'LineWidth'},{2;2;2})
             xlabel('days'); ylabel('STB');
             title(['Water production: producer ', num2str(i)])
             legend('Reference','Initial model','Calibrated model','Location','SouthEast')     
              hl = line([T(end) T(end)],get(hs,'YLim'),'Color', 150*[1 1 1]/255,'LineStyle','--','HandleVisibility','off');
            set(hl,'LineWidth',lineSize)
  end
  
  figure
  
  for i = 1:4
     hs =  subplot(2,2,i)
           h5 = plot( hs ,tt,-qo_ref(:,i)*f_units,'.k',tt,-qo_0(:,i)*f_units,tt,-qo_pred(:,i)*f_units)
             set(h5,{'LineWidth'},{2;2;2})
             xlabel('days'); ylabel('STB');
             title(['Oil production: producer ', num2str(i)])
             legend('Reference','Initial model','Calibrated model','Location','NorthEast')     
             hl = line([T(end) T(end)],get(hs,'YLim'),'Color', 150*[1 1 1]/255,'LineStyle','--','HandleVisibility','off');
             set(hl,'LineWidth',lineSize)
  end
 
  
 
 
  
  %% BHP injectors
  
 for i = 1:numel(schedule.step.val)
     for p = 1:8
        bhp_ref(i,p) = wellSols_ref{i}(p).bhp;
        bhp_0(i,p) = wellSols_0{i}(p).bhp;
        bhp_pred(i,p) = wellSols_opt{i}(p).bhp;
     end
 end
 tt = cumsum(dt)/day;
 f_units = 1/barsa;

  
  II = [1:8]
  figure
   for i = 1:8
            hs = subplot(2,4,II(i))
            h7 = plot(hs , tt,bhp_ref(:,i)*f_units,'.k',tt, bhp_0(:,i)*f_units,tt, bhp_pred(:,i)*f_units)
             set(h7,{'LineWidth'},{2;2;2})
             xlabel('days'); ylabel('bar');
             title(['bhp: injector ', num2str(i)])
             legend('Reference','Initial model','Calibrated model','Location','SouthEast')  
             hl = line([T(end) T(end)],get(hs,'YLim'),'Color', 150*[1 1 1]/255,'LineStyle','--','HandleVisibility','off');
             set(hl,'LineWidth',lineSize)
   end
   
  