function [val, der, wellSols, states, extra] = evalObjectiveAndSens(u, obj, state0, model, schedule, scaling,varargin)
% Objective (and gradient) evaluation function based on input control vector u
% u(1:n) : dz
% u(n+1) : rhomult
% u(n+2) : permmult
% u(n+3)  : poromult
opt=struct('only_sens',false,...
           'wellSols',[],...
           'states',[]);
opt = merge_options(opt,varargin{:});

% where minimizing -> set objective to minus objective:

minu = min(u);
maxu = max(u);
if or(minu < -eps , maxu > 1+eps)
    warning('Controls are expected to lie in [0 1]\n')
end

boxLims = scaling.boxLims;
if isfield(scaling, 'obj')
    objScaling = scaling.obj;
else
    objScaling = 1;
end

% update model
% dz, rhofac, permfac, porofac

n = model.G.cells.num;
% 1.

% multipliers
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
us = u.*(umax-umin)+umin;
q_prev = arrayfun(@(x)x.W.val, schedule.control);

for i=1:numel(schedule.control)
    schedule.control(i).dz       = us(1:n);
    schedule.control(i).W.val    = us(n+1)*schedule.control(i).W.val/mean(q_prev);
    schedule.control(i).rhofac   = us(n+2);
    schedule.control(i).permfac  = us(n+3);
    schedule.control(i).porofac  = us(n+4);
end
if(nargout==5)
  extra.schedule=schedule;
  extra.state0=state0;
end    
% run simulation:
if(opt.only_sens)
  assert(~isempty(opt.wellSols))
  wellSols=opt.wellSols;
  states=opt.states;    
else
[wellSols, states, sim_report] = simulateScheduleAD(state0, model, schedule);
end

states = addHeightData(states, model.G, model.fluid);
% compute objective:
vals = obj(wellSols, states, schedule);
vals=horzcat(vals{:})';
val  = - sum(vals,1)/objScaling;
%val=nan(1,numel(vals));
%for i=1:numel(vals)
%    val(i)  = - sum(cell2mat(vals{i}))/objScaling;
%end
% run adjoint:
if nargout > 1
    objh = @(tstep)obj(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
    g    = computeGradientAdjointAD(state0, states, model, schedule, objh, 'ControlVariables', {'scell','well','mult'});    
    %g    = cell2mat(g);
    % need to take special care of rate-multiplier
    for i=1:numel(q_prev)
       g{2,i}= g{2,i}.*q_prev(i)/mean(q_prev);
    end
    %g(n+1,:) = g(n+1,:).*q_prev/mean(q_prev);
    if(nargout==5)
        extra.der=g;
    end
    gsum=vertcat(g{:,1});
    for i=2:size(g,2)
        gsum=gsum+vertcat(g{:,i});
    end
    g=gsum;
    %g    = sum(g,2);
    % scale gradient:
    dBox   = boxLims(:,2) - boxLims(:,1);
    der  = - bsxfun(@times,g,(dBox/objScaling));
    if(nargout==5)
       derold=extra.der;
       extra.der={};
       for i=1:size(extra.der,2) 
        extra.der{i}=bsxfun(@times,vertcat(derold{:,i}),-dBox/objScaling);
       end
    end
    
end
end

% function grd = scaleGradient(grd, schedule, boxLims, objScaling)
% dBox   = boxLims(:,2) - boxLims(:,1);
% for k = 1:numel(schedule.control)
%     grd{k} = (dBox/objScaling).*grd{k};
% end
% end
    
