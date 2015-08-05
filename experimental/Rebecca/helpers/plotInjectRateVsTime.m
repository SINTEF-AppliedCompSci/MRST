function [ hfig, hax ] = plotInjectRateVsTime(schedule,inj_year,rhoCref)

    figure; set(gcf,'Position',[1000 1000 1500 400])
    
    subplot(1,4,1)
    timeSinceInjLast        = 0;
    [timeSinceInj, rateNow] = deal(zeros(1,numel(schedule.step.val)));
    for i = 1:numel(schedule.step.val)
        timeSinceInj(i)  = schedule.step.val(i)/365/24/60/60 + timeSinceInjLast; % years
        rateNow(i)       = schedule.control(schedule.step.control(i)).W.val; % m^3/s
        timeSinceInjLast = timeSinceInj(i);
    end
    plot(timeSinceInj, rateNow, 'o--k')
    xlabel('time, yr','FontSize',14)
    ylabel('injection rate, m^3/s','FontSize',14)
    hold off
    ax = gca;
    ax.XTickLabel = ax.XTick + inj_year(1)-1;

    subplot(1,4,2)
    plot(timeSinceInj, rateNow*(rhoCref/1e9*365*24*60*60), 'o--k') % Mt/yr
    xlabel('time, yr','FontSize',14')
    ylabel('injection rate, Mt/year','FontSize',14)
    hold off
    ax = gca;
    ax.XTickLabel = ax.XTick + inj_year(1)-1;

    subplot(1,4,3) %(compare this plot against Cavanagh 2013, fig 3)
    massNowLast             = 0;
    [timeStepNow, massNow]  = deal(zeros(1,numel(schedule.step.val)));
    for i = 1:numel(schedule.step.val)
        timeStepNow(i)  = schedule.step.val(i)/(365*24*60*60); % yr
        massNow(i)      = rateNow(i)*(rhoCref/1e9*365*24*60*60)*timeStepNow(i) + massNowLast; % Mt
        massNowLast     = massNow(i);
    end
    plot(timeSinceInj, massNow, 'o--k')
    xlabel('time, yr','FontSize',14)
    ylabel('Accumulated mass injected, Mt','FontSize',14)
    hold off
    ax = gca;
    ax.XTickLabel = ax.XTick + inj_year(1)-1;
    
    
    subplot(1,4,4) % compare this plot with given Sleipner injection rate data (mass, kg).
    massNow_kg = massNow.*1e9; % kg
    plot(timeSinceInj, massNow_kg, 'o--k')
    xlabel('time, yr','FontSize',14)
    ylabel('Accumulated mass injected, kg','FontSize',14)
    hold off
    ax = gca;
    ax.XTickLabel = ax.XTick + inj_year(1)-1;
    
    hfig = gcf;
    hax  = gca;

end

