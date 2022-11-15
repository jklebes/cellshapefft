function varargout = Main_Gui(varargin)
% MAIN_GUI MATLAB code for Main_Gui.fig
%      MAIN_GUI, by itself, creates a new MAIN_GUI or raises the existing
%      singleton*.
%
%      H = MAIN_GUI returns the handle to a new MAIN_GUI or the handle to
%      the existing singleton*.
%
%      MAIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GUI.M with the given input arguments.
%
%      MAIN_GUI('Property','Value',...) creates a new MAIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_Gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_Gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main_Gui

% Last Modified by GUIDE v2.5 30-May-2018 18:29:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_Gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_Gui_OutputFcn, ...
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


% --- Executes just before Main_Gui is made visible.
function Main_Gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main_Gui (see VARARGIN)

% Choose default command line output for Main_Gui
handles.output = hObject;
handles.Param = struct('siz',[],'tleng',[],'pas2',[],'name',[],'fres',[]...
    ,'cut',[],'propor',[],'pathin',[],'pathout',[],'timestep',[]...
    ,'sigma',[],'strel',[],'register',[],'regsize',[],'nbpoints',[],'sizeview',[],'contour',[],'rec',[],'tsart',[]);
handles.Param.scale = 150;                    % scaling for cell deformation 
handles.Param.fres = 'results';
handles.Param.contour=[];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Main_Gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_Gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse1.
function browse1_Callback(hObject, eventdata, handles)
% hObject    handle to browse1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of browse1
lala =  uigetdir();
handles.Param.pathout = [lala,'/', handles.Param.fres ,'_averageon_', num2str(handles.Param.timestep) '/']; % path to result folder
if ~exist(handles.Param.pathout,'dir')             % creation of result folder
mkdir(handles.Param.pathout)
end
handles.Param.pathin = '';
%set(handles.pathmain,'String',[handles.Param.pathin,'/']);

guidata(hObject, handles);




% --- Executes on button press in browse2.
function browse2_Callback(hObject, eventdata, handles)
% hObject    handle to browse2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of browse2
handles.Param.name = uigetdir2();
set(handles.pathim,'String','Done !');
guidata(hObject, handles);



function pathmain_Callback(hObject, eventdata, handles)
% hObject    handle to pathmain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pathmain as text
%        str2double(get(hObject,'String')) returns contents of pathmain as a double
handles.Param.pathin = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pathmain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathmain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pathim_Callback(hObject, eventdata, handles)
% hObject    handle to pathim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pathim as text
%    handles.Param.pathin = get(hObject,'String');    str2double(get(hObject,'String')) returns contents of pathim as a double
handles.Param.name = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pathim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pathim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushOK.
function pushOK_Callback(hObject, eventdata, handles)
% hObject    handle to pushOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Param.pathin = [handles.Param.pathin,'/'];
for kk = 1:numel(handles.Param.name)
    handles.Param.name{kk} = [handles.Param.name{kk},'/'];
end
Parameters = handles.Param;
%save([handles.Param.pathout,'Param.mat'],'Parameters');
GUI_deformation_ft(Parameters)
close(Main_Gui());




function nb_Callback(hObject, eventdata, handles)
% hObject    handle to nb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nb as text
%        str2double(get(hObject,'String')) returns contents of nb as a double
handles.Param.tleng = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbave_Callback(hObject, eventdata, handles)
% hObject    handle to nbave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbave as text
%        str2double(get(hObject,'String')) returns contents of nbave as a double
handles.Param.timestep = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nbave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function recouvrement_Callback(hObject, eventdata, handles)
% hObject    handle to recouvrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recouvrement as text
%        str2double(get(hObject,'String')) returns contents of recouvrement as a double
handles.Param.rec = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function recouvrement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recouvrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTSART_Callback(hObject, eventdata, handles)
% hObject    handle to editTSART (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTSART as text
%        str2double(get(hObject,'String')) returns contents of editTSART as a double
handles.Param.tsart = str2double(get(hObject,'String'));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function editTSART_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTSART (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editscale_Callback(hObject, eventdata, handles)
% hObject    handle to editscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editscale as text
%        str2double(get(hObject,'String')) returns contents of editscale as a double
handles.Param.scale = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
