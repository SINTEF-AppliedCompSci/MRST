function ShowSlider(model)
%
% DESCRIPTION: adds the time slides to the plots created with the Visualize
% module
%
% SYNOPSIS:
%   ShowSlider(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - state: states of the simulation in MRST format
%   - grid
%   - dynamic: information saved during the simulation
%   - experiment: saturation functions used for forward modeling
%   - simulation: time stepping and grid cells information
%
% RETURNS:
%   time slider added to the simulation prediction plots
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
    static = model.static;    
    fig = static.fig;
    figTitle = fig.title;
    if isfield(fig, 'style')
        figStyle    = fig.style;
    else
        figStyle = 'docked';
    end
    handles.f = findobj('type','figure','name',figTitle);
    handles.ax = findobj(handles.f,'type','axes');
    dynamic = model.dynamic;
    params  = dynamic.params;
    sliderValue = params.counter;    
    handles.model = model;
    axs = handles.ax;
    title1 = fig.subtitle{1};
    if(strcmp(axs(1).Title.String,strcat(title1))), sliderAx = axs(1); end
    if(strcmp(axs(2).Title.String,strcat(title1))), sliderAx = axs(2); end
    if(strcmp(axs(3).Title.String,strcat(title1))), sliderAx = axs(3); end
    if(strcmp(axs(4).Title.String,strcat(title1))), sliderAx = axs(4); end
    if(strcmp(axs(5).Title.String,strcat(title1))), sliderAx = axs(5); end
    if(strcmp(axs(6).Title.String,strcat(title1))), sliderAx = axs(6); end
    figPos = get(sliderAx, 'position');
    if(strcmp(figStyle,'normal'))
        sliderW = figPos(3)/1.3;
        sliderH = figPos(end)/7;
        sliderOff = [0.11 -0.06];
    else
        sliderW = figPos(3)/1.4;
        sliderH = figPos(end)/7;
        sliderOff = [0.09 -0.08];
    end
    sliderPos = [figPos(1:2)+sliderOff sliderW sliderH];
     
%     handles.editbox = uicontrol('style','edit');
    handles.slider  = uicontrol('style','slider');
    handles.textbox = uicontrol('style','text');
    handles.popup   = uicontrol('style','popupmenu');
        
    % set slider properties
    set(handles.slider,'units',get(sliderAx,'units'))    
    set(handles.slider,'position',sliderPos);
    set(handles.slider,'min',1,'max',sliderValue);
    set(handles.slider,'value', sliderValue);
    set(handles.slider,'sliderstep', [1/(sliderValue-1) 1/(sliderValue-1)]);
    set(handles.slider,'callBack',{@SliderCallback,handles});
    
    % set box properties
%     boxW = sliderH;
%     boxH = sliderH;
%     boxOff = [0.005 0 0];      
%     boxPos = [sliderPos(1)+sliderW+boxOff(1) sliderPos(2)+boxOff(2) boxW boxH];              
%     set(handles.editbox,'units',get(sliderAx,'units'))
%     set(handles.editbox,'position', boxPos)
%     set(handles.editbox,'string', sliderValue)
%     set(handles.editbox,'callback',{@EditboxCallback})
    
    % set label text properties
    textW = 2 * sliderH;
    textH = sliderH;
    textOff = [0.005 0.03 0];      
    textPos = [sliderPos(1)+sliderW+textOff(1) sliderPos(2)+textOff(2) textW textH];              
    set(handles.textbox,'units',get(sliderAx,'units'))
    set(handles.textbox,'position', textPos)
    % display unit
    displayUnits = DisplayUnits(model);    
    displayTime  = displayUnits.displayTime;
    set(handles.textbox,'string',strcat('Time (',displayTime,')'),'fontsize',8.5);
    time = params.cumScheduleSteps;
    timeFac = Convert(displayTime);
    time = time ./ timeFac;
    popupValues = round(time,3);
    
    % set popup properties
    popupW = 2 * sliderH;
    popupH = sliderH;
    popupOff = [0.005 0 0];      
    popupPos = [sliderPos(1)+sliderW+popupOff(1) sliderPos(2)+popupOff(2) popupW popupH];              
    set(handles.popup,'units',get(sliderAx,'units'))
    set(handles.popup,'position', popupPos)
    set(handles.popup,'string', popupValues)
    set(handles.popup,'callback',{@PopupCallback})
     
    % store s struct in the figure
    guidata(handles.f, handles);
end

function SliderCallback(hObject, eventdata, handles)
    sliderValue = round(get(hObject,'value'));
%     set(handles.editbox,'string', num2str(sliderValue));    
    model = handles.model;
    dynamic = model.dynamic;
    params  = dynamic.params;
    set(handles.popup,'value',sliderValue);
    params.fromIdx = 1;
    params.toIdx = sliderValue; 
    dynamic.params = params;
    model.dynamic  = dynamic;
    UpdatePlots(model);
end

function EditboxCallback(hObject, eventdata)
    handles = guidata(hObject);
    sliderValue = str2double(get(hObject,'string'));    
    model = handles.model;
    dynamic = model.dynamic;
    params  = dynamic.params;
    if(sliderValue > params.counter)
        sliderValue = params.counter;
        set(handles.editbox,'string', num2str(sliderValue));
    end    
    set(handles.slider,'value',sliderValue);
    set(handles.popup,'string',params.cumScheduleSteps(sliderValue));
    params.fromIdx = 1;
    params.toIdx = sliderValue;
    dynamic.params = params;
    model.dynamic  = dynamic;
    UpdatePlots(model);
end

function PopupCallback(hObject, eventdata)
    handles = guidata(hObject); 
    idx = get(hObject,'Value');   
    model = handles.model;
    dynamic = model.dynamic;
    params  = dynamic.params;  
    set(handles.slider,'value',idx);
    params.fromIdx = 1;
    params.toIdx = idx;
    dynamic.params = params;
    model.dynamic  = dynamic;
    UpdatePlots(model);
end
