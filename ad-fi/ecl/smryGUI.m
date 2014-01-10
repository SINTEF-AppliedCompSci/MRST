function varargout = smryGUI(varargin)
% SMRYGUI M-file for smryGUI.fig
%      SMRYGUI, by itself, creates a new SMRYGUI or raises the existing
%      singleton*.
%
%      H = SMRYGUI returns the handle to a new SMRYGUI or the handle to
%      the existing singleton*.
%
%      SMRYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMRYGUI.M with the given input arguments.
%
%      SMRYGUI('Property','Value',...) creates a new SMRYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before smryGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to smryGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help smryGUI

% Last Modified by GUIDE v2.5 29-Oct-2012 17:38:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @smryGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @smryGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before smryGUI is made visible.
function smryGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to smryGUI (see VARARGIN)

if nargin < 4
    if isempty(evalin('base', 'who(''smry'')'))
        error('No smry-variable');
    end
    handles.smry = evalin('base', 'smry');
else
    handles.smry = varargin{1};
end

global smry;
smry = handles.smry;


nm   = handles.smry.WGNAMES{1};
kws  = handles.smry.getKws(nm);
handles.extFig = false;
handles.cur_nm   = nm;
handles.cur_kw  = kws{1};
handles.cur_ix  = find( handles.smry.getInx(handles.cur_nm, handles.cur_kw) );
tn = handles.smry.getNms('TIME');
handles.ms = handles.smry.get(tn(1),'TIME',:);
handles.plotFormat = '-b';
handles.plotIx     = (1:numel(handles.cur_ix))';
makeplot(handles)

% Choose default command line output for smryGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes smryGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = smryGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
str = get(hObject, 'String');
val = get(hObject, 'Value');
cur_nm =  str{val};
handles.cur_nm = cur_nm;
kws = handles.smry.getKws(cur_nm);
[tf, loc] = ismember(handles.cur_kw, kws);
if loc == 0, loc = 1;end
set(handles.popupmenu2, 'String', kws, 'Value', loc);
if ~ismember(handles.cur_kw, handles.smry.getKws(cur_nm))
    handles.cur_kw = kws{1};
end
handles.cur_ix  = find( handles.smry.getInx(handles.cur_nm, handles.cur_kw) );
ixStr = vertcat({'all'}, cellstr(num2str( (1:numel(handles.cur_ix))' )));
set(handles.popupmenu4, 'String', ixStr, 'Value', 1);
handles.plotIx = 0;
guidata(hObject,handles)
makeplot(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global smry
set(hObject, 'String', smry.WGNAMES);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
str = get(hObject, 'String');
val = get(hObject, 'Value');
cur_kw =  str{val};
handles.cur_kw = cur_kw;

handles.cur_ix  = find( handles.smry.getInx(handles.cur_nm, handles.cur_kw) );
ixStr = vertcat({'all'}, cellstr(num2str( (1:numel(handles.cur_ix))' )));
set(handles.popupmenu4, 'String', ixStr, 'Value', 1);
handles.plotIx = 0;

guidata(hObject,handles)
makeplot(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global smry
set(hObject, 'String', smry.WGNAMES{1});
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
handles.extFig = get(hObject,'Value');
guidata(hObject,handles);
makeplot(handles)



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function makeplot(handles)
if ismember(handles.cur_nm, handles.smry.getNms(handles.cur_kw))
    x = handles.ms;
    y = handles.smry.get(handles.cur_nm, handles.cur_kw,:);
    if handles.plotIx > 0
        y = y(handles.plotIx, :);
    end
    if handles.extFig
        figure(1);
    end
    plot(x,y,handles.plotFormat);
else
    plot(-1);
    axis([0 handles.ms(end) -1 1]);
end
title(handles.cur_nm)
ylabel(handles.cur_kw)
xlabel('time [days]');



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.plotFormat = get(hObject, 'String');
guidata(hObject,handles);
makeplot(handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
%str = get(hObject, 'String');
val = get(hObject, 'Value');
handles.plotIx = val-1;
guidata(hObject,handles)
makeplot(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
set(hObject, 'String', 'all');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
