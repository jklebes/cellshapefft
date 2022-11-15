function varargout = GUI_deformation_ft(varargin)
% GUI_deformation_ft MATLAB code for GUI_deformation_ft.fig
%      GUI_deformation_ft, by itself, creates a new GUI_deformation_ft or raises the existing
%      singleton*.
%
%      H = GUI_deformation_ft returns the handle to a new GUI_deformation_ft or the handle to
%      the existing singleton*.
%
%      GUI_deformation_ft('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_deformation_ft.M with the given input arguments.
%
%      GUI_deformation_ft('Property','Value',...) creates a new GUI_deformation_ft or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_deformation_ft_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_deformation_ft_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_deformation_ft

% Last Modified by GUIDE v2.5 18-Feb-2019 17:32:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_deformation_ft_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_deformation_ft_OutputFcn, ...
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
% 
% function cleanednames = cleanfolder(folder)
%     int1 = dir(folder);
%     int1([int1.isdir]) = [];
%     A = strfind({int1.name},'.DS_Store');
%     cleanednames = int1(cellfun('isempty',A));


% --- Executes just before GUI_deformation_ft is made visible.
function GUI_deformation_ft_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_deformation_ft (see VARARGIN)

% Choose default command line output for GUI_deformation_ft
handles.output = hObject;

handles.Param=varargin{1};

handles.Anamesi = cleanfolder([handles.Param.name{1}]);

handles.image1 = 0;
handles.image2 =0;
handles.image3 =0;
handles.Param.regsize = 0;
handles.Param.cut = 0;
handles.Param.register = 0;
handles.Param.propor = 0.1;
handles.Param.sigma = 0.1;
handles.Param.strel = 4;
handles.Param.pas2 = 128;
handles.Param.zoom = 6;

handles.rec = (1-handles.Param.rec)*handles.Param.pas2;
handles.Param.nbpoints = 30;
handles.Param.sizeview = 0;

handles = calculregion(handles);
handles = affichim(handles);
handles = fftcalculus(handles);
handles = gaussfilt(handles);
handles = cutfun(handles);
handles = affichcutfft(handles);
handles = proporfft(handles);
handles = affichproporfft(handles);
handles = filter(handles);
handles = trace_defo(handles);
set(get(handles.filt1,'children'),'visible','off')
set(get(handles.filt2,'children'),'visible','off')
set(get(handles.filt3,'children'),'visible','off')

KB=handles.Param.pas2^2*9.53674e-7*8;  
KBwhole  = KB*handles.regl;  
set(handles.sizesub,'String',sprintf(['One image: ' num2str(floor(KBwhole)) ' MB']));
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_deformation_ft_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;

% --- Executes on button press in shuffle.
function shuffle_Callback(hObject, eventdata, handles)
% hObject    handle to shuffle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shuffle
shuffleval = get(hObject,'Value');
    if shuffleval == 1
        set(handles.shuffle,'Value',0);
        handles = calculregion(handles);
        handles = affichim(handles);
        handles = fftcalculus(handles);
        handles = gaussfilt(handles);
        handles = cutfun(handles);
        handles = affichcutfft(handles);
        handles = proporfft(handles);
        handles = affichproporfft(handles);
        handles = trace_defo(handles);
        
        if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

        
        if handles.adv == 1
            handles = filter(handles);
        else
            set(get(handles.filt1,'children'),'visible','off')
            set(get(handles.filt2,'children'),'visible','off')
            set(get(handles.filt3,'children'),'visible','off')

        end
    else
    end
guidata(hObject, handles);



% --- Executes on button press in choose.
function choose_Callback(hObject, eventdata, handles)
% hObject    handle to choose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of choose
     a = imread([handles.Param.name,  handles.Anamesi(handles.Param.tsart).name]);
     a = im2double(a);

if get(hObject,'Value')
    F = figure(2);
    imshow(a)
    hold on 
    for win = 1:2:handles.regl
       x = handles.Posi(win,1)+handles.Param.pas2/2;
       y = handles.Posi(win,2)+handles.Param.pas2/2;
       text(x,y,num2str(win),'Color','yellow','FontSize',5);
    
    end
    Numbs = choose();
    waitfor(Numbs.OK,'Value',1)
    handles.image1 = str2double(Numbs.NB1.String);   
    handles.image2 = str2double(Numbs.NB2.String);
    handles.image3 = str2double(Numbs.NB3.String);
    close(Numbs.figure1);
%        handles.image1 = input('What is the number of the first subimage? ');
%        handles.image2 = input('What is the number of the second subimage? ');
%        handles.image3 = input('What is the number of the third subimage? ');
     close(F); 
        handles = calculregion(handles);
        handles = affichim(handles);
        handles = fftcalculus(handles);
        handles = gaussfilt(handles);
        handles = cutfun(handles);
        handles = affichcutfft(handles);
        handles = proporfft(handles);
        handles = affichproporfft(handles);
        handles = trace_defo(handles);
     if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

        if handles.adv == 1
            handles = filter(handles);
        else
            set(get(handles.filt1,'children'),'visible','off')
            set(get(handles.filt2,'children'),'visible','off')
            set(get(handles.filt3,'children'),'visible','off')

        end
else
    handles.image1 = 0;
    handles.image2 =0;
    handles.image3 =0;
        handles = calculregion(handles);
        handles = affichim(handles);
        handles = fftcalculus(handles);
        handles = gaussfilt(handles);
        handles = cutfun(handles);
        handles = affichcutfft(handles);
        handles = proporfft(handles);
        handles = affichproporfft(handles);
        handles = trace_defo(handles);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

        
        if handles.adv == 1
            handles = filter(handles);
        else
            set(get(handles.filt1,'children'),'visible','off')
            set(get(handles.filt2,'children'),'visible','off')
            set(get(handles.filt3,'children'),'visible','off')

        end
end
guidata(hObject, handles);





% --- Executes on slider movement.
function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%newval = floor(get(hObject,'Value'));
    newval = get(hObject,'Value');
 %   handles.sigma = newval;
    handles = fftcalculus(handles);
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);
    set(handles.editsigma,'String',newval);
    handles.Param.sigma = newval;
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

   guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
set(hObject,'Min',0,'Max',2,'SliderStep',[0.1,1]);
%set(hObject,'Min',0,'Max',100,'SliderStep',[1,10],'Value',1);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function propor_Callback(hObject, eventdata, handles)
% hObject    handle to propor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    handles.Param.propor = get(hObject,'Value');
    handles.Param.propor = handles.Param.propor*10^-3;
    handles = proporfft(handles);
    handles = affichproporfft(handles);
    handles = trace_defo(handles);
    set(handles.editpropor,'String',handles.Param.propor);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
    if get(handles.adv,'Value')
    handles = filter(handles);
    else
    end
   
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function propor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.

set(hObject,'Min',0,'Max',200,'SliderStep',[1/50,10/50],'Value',10);

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject,handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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
 guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function strel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
set(hObject,'Min',5,'Max',10,'SliderStep',[0.2,2],'Value',5,'Visible','off');

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject,handles);



% --- Executes on slider movement.
function strel_Callback(hObject, eventdata, handles)
% hObject    handle to strel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    newval = floor(get(hObject,'Value'));
    if newval > 4 && newval ~= handles.Param.strel
        handles.Param.strel = newval;
        handles = proporfft(handles);
        handles = affichproporfft(handles);
        handles = trace_defo(handles);
        set(handles.editstrel,'String',handles.Param.strel);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
   handles = filter(handles);

    else
    end
guidata(hObject,handles);

% 
% % --- Executes on button press in pushsigma.
% function pushsigma_Callback(hObject, eventdata, handles)
% % hObject    handle to pushsigma (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.pushsigma,'String',handles.sigma);
% handles.Param.sigma = handles.sigma;
% %save([handles.Param.pathout,'Param.mat'],'handles.Param');
% guidata(hObject,handles);


% % --- Executes on button press in pushstrel.
% function pushstrel_Callback(hObject, eventdata, handles)
% % hObject    handle to pushstrel (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.pushstrel,'String',handles.Param.strel);
% %save([handles.Param.pathout,'Param.mat'],'handles.Param');
% guidata(hObject,handles);


% --- Executes on slider movement.
function window_Callback(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

nl = round(get(hObject,'Value'));
% n = log(handles.Param.pas2)/log(2);
% if nl == n
% else
   % set(hObject,'Value',nl);
%    handles.Param.pas2 = 2^nl;
 
%handles.rec = 1/2*handles.Param.pas2;
 
    handles.Param.pas2 = nl;
    if floor(handles.Param.pas2/2)==handles.Param.pas2/2
        handles.rec = handles.Param.pas2/2;
     else
         handles.rec = (handles.Param.pas2-1)/2;
    end

    handles = calculregion(handles);
    handles = affichim(handles);
    KB=handles.Param.pas2^2/16*9.53674e-7*8;  
    KBwhole = KB*handles.regl;
    set(handles.sizesub,'String',sprintf(['One image: ' num2str(floor(KBwhole)) ' MB']))
    handles = affichim(handles);
    handles = fftcalculus(handles);
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);
    handles = proporfft(handles);
    handles = affichproporfft(handles);
    handles = trace_defo(handles);
    set(handles.enterwindow,'String',handles.Param.pas2);
    
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

    if handles.adv == 1
        handles = filter(handles);
    else
        set(get(handles.filt1,'children'),'visible','off')
        set(get(handles.filt2,'children'),'visible','off')
        set(get(handles.filt3,'children'),'visible','off')

    end
%end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
set(hObject,'Min',20,'Max',1000,'Value',100,'SliderStep',[2/980,20/980]);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject,handles);

% 
% % --- Executes on button press in pushwindow.
% function pushwindow_Callback(hObject, eventdata, handles)
% % hObject    handle to pushwindow (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.pushwindow,'String',handles.Param.pas2);
% guidata(hObject,handles);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


   
% % --- Executes on button press in pushcut.
% function pushcut_Callback(hObject, eventdata, handles)
% % hObject    handle to pushcut (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.pushcut,'String',get(handles.cut,'Value'));
% handles.Param.cut = get(handles.cut,'Value');
%    
% %save([handles.Param.pathout,'Param.mat'],'handles.Param');
% guidata(hObject,handles);
% %guidata(handles,hObject);
% 
% % --- Executes on button press in pushpropor.
% function pushpropor_Callback(hObject, eventdata, handles)
% % hObject    handle to pushpropor (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.pushpropor,'String',handles.propor);
% handles.Param.propor = handles.propor;
% %save([handles.Param.pathout,'Param.mat'],'handles.Param');
% guidata(hObject,handles);


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Val = get(hObject,'Value');
if Val == 1
   % Parameters = struct('propor',handles.propor,'cut',handles.cut,'strel',handles.strel,'sigma',handles.sigma,'pas2',handles.Param.pas2,'registering',handles.register);
    Parameters=handles.Param;
    save([handles.Param.pathout,'Param.mat'],'Parameters');
    
   savefig(handles.figure1,[handles.Param.pathout,'Parameters_fig.fig'])     % save image fig

    hgexport(handles.figure1,[handles.Param.pathout,'Parameters_fig.png'],hgexport('factorystyle'),'Format','png');
    
    close(handles.figure1);

    full_analysis(handles.Param);
else
end



% Hint: get(hObject,'Value') returns toggle state of OK



% --- Executes on button press in sizes.
function sizes_Callback(hObject, eventdata, handles)
% hObject    handle to sizes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sizes
handles.Param.sizeview = get(hObject,'Value');


if get(hObject,'Value')  && get(handles.anisotrpy_ellipse,'Value')==0
    handles = ellipses_def(handles);
    handles = affichsize(handles);
end
    if get(hObject,'Value')  && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
    if get(hObject,'Value') ==0 && get(handles.anisotrpy_ellipse,'Value')==0 
    handles = affichcutfft(handles);
    handles = trace_defo(handles);
    end
guidata(hObject,handles);




% --- Executes on button press in regsize.
function regsize_Callback(hObject, eventdata, handles)
% hObject    handle to regsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of regsize
handles.Param.regsize = get(hObject,'Value');

guidata(hObject, handles);




% --- Executes on slider movement.
function slidersize_Callback(hObject, eventdata, handles)
% hObject    handle to slidersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.Param.nbpoints = get(hObject,'Value');
set(handles.editsize,'String',handles.Param.nbpoints);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
    if get(handles.sizes,'Value')==0 && get(handles.anisotrpy_ellipse,'Value')==0 
    handles = affichcutfft(handles);
    handles = trace_defo(handles);
    end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slidersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
set(hObject,'Min',10,'Max',200,'SliderStep',[10/190,20/190],'Value',10,'Visible','off');

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on button press in adv.
function adv_Callback(hObject, eventdata, handles)
% hObject    handle to adv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of adv
if get(hObject,'Value')
    set(handles.editstrel,'Visible','On')
    set(handles.regsize,'Visible','On')
    set(handles.textstrel,'Visible','On')
    set(handles.strel,'Visible','On')
    set(handles.text11,'Visible','On')
    set(get(handles.filt1,'children'),'visible','on')
    set(get(handles.filt2,'children'),'visible','on')
    set(get(handles.filt3,'children'),'visible','on')
    set(handles.text18,'visible','on')
    set(handles.sizes,'visible','on')
    set(handles.slidersize,'visible','on')
     set(handles.editsize,'visible','on')
     set(handles.anisotrpy_ellipse,'Visible','On')
    set(handles.nameani_ellipse,'Visible','On')

     handles = filter(handles);

else
    set(handles.editstrel,'Visible','Off')
      set(handles.regsize,'Visible','Off')
    set(handles.textstrel,'Visible','Off')
    set(handles.strel,'Visible','Off')
     set(handles.text18,'visible','Off')
       set(handles.anisotrpy_ellipse,'Visible','Off')
    set(handles.nameani_ellipse,'Visible','Off')

    set(handles.text11,'Visible','Off')
    set(get(handles.filt1,'children'),'visible','off')
    set(get(handles.filt2,'children'),'visible','off')
    set(get(handles.filt3,'children'),'visible','off')
    set(handles.sizesstring,'visible','off')
    set(handles.sizes,'visible','off')   
    set(handles.slidersize,'visible','off')
 set(handles.editsize,'visible','off')
   
end
guidata(hObject,handles);



function handles=calculregion(handles)
 
[handles.Posi,handles.regl] = Position_gui(handles.Param);

function handles=fftcalculus(handles)
% asubs contains the subimages we will use 
     FTs = avfft(handles.Param.name,handles.wins,handles.Param,handles.Posi);
     handles.fftsub1 = reshape(FTs(1,:,:),handles.Param.pas2,handles.Param.pas2);
    handles.fftsub2 = reshape(FTs(2,:,:),handles.Param.pas2,handles.Param.pas2);
    handles.fftsub3 = reshape(FTs(3,:,:),handles.Param.pas2,handles.Param.pas2);

     
function handles = gaussfilt(handles)
    filtre = gaussianFilter(ceil(2*handles.Param.sigma),handles.Param.sigma);
    handles.spsub1 = filter2(filtre,handles.fftsub1);
    handles.spsub2 = filter2(filtre,handles.fftsub2);
    handles.spsub3 = filter2(filtre,handles.fftsub3);
   
        

% --- Executes on slider movement.
function cut_Callback(hObject, eventdata, handles)
% hObject    handle to cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newval = floor(get(hObject,'Value'));
if  newval ~= handles.Param.cut
    set(hObject,'Value',newval);
    handles.Param.cut=newval;
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);
    set(handles.editcut,'String',handles.Param.cut);
    

    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
        
else
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
set(hObject,'Min',0,'Max',30,'SliderStep',[0.1,1],'Value',1);

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject,handles);


% --- Executes on button press in register.
function register_Callback(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of register
handles.Param.register = get(hObject,'Value');
guidata(hObject, handles);




function sizesub_Callback(hObject, eventdata, handles)
% hObject    handle to sizesub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sizesub as text
%        str2double(get(hObject,'String')) returns contents of sizesub as a double


% --- Executes during object creation, after setting all properties.
function sizesub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sizesub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = cutfun(handles)
    handles.fftsubcut1 = handles.spsub1;
    handles.fftsubcut2 = handles.spsub2;
    handles.fftsubcut3 = handles.spsub3;
 ic = handles.Param.pas2*1/2+1;
 jc = handles.Param.pas2*1/2+1;
 [I,J] = ndgrid(1:handles.Param.pas2,1:handles.Param.pas2);
 M = double((I-ic).^2+(J-jc).^2<=(handles.Param.cut)^2); 
 handles.fftsubcut1(M>0) = 0;
 handles.fftsubcut2(M>0) = 0;
 handles.fftsubcut3(M>0) = 0;
 
 
 
function handles=affichcutfft(handles)

   imagesc(handles.fftsubcut1,'Parent',handles.fft1cut)
   axes(handles.fft1cut);
   zoom(handles.Param.zoom)
   set(handles.fft1cut,'Visible','off')

   imagesc(handles.fftsubcut2,'Parent',handles.fft2cut)
   axes(handles.fft2cut);
   zoom(handles.Param.zoom)
   set(handles.fft2cut,'Visible','off')

   imagesc(handles.fftsubcut3,'Parent',handles.fft3cut)
   axes(handles.fft3cut);
   zoom(handles.Param.zoom)
   set(handles.fft3cut,'Visible','off')


function handles=affichim(handles)
     
    a = choose_image(handles.Param.name{1},handles.Anamesi,handles.Param.tsart);

     if handles.image1 ~=0
         win1 = handles.image1;
         win2 = handles.image2;
         win3 = handles.image3;
     elseif handles.regl == 1
         win1 = 1;
         win2 = 1;
         win3 = 1;         
     else
       win1 = round(rand(1)*handles.regl);
       win2 = round(rand(1)*handles.regl);
       win3 = round(rand(1)*handles.regl);
         
     end
     
     set(handles.num1,'String',['Number: ' num2str(win1)]);
     set(handles.num2,'String',['Number: ' num2str(win2)]);
     set(handles.num3,'String',['Number: ' num2str(win3)]);
     
     x = handles.Posi(win1,1);
     y = handles.Posi(win1,2);
     
     handles.sub1=a(y:y+handles.Param.pas2-1,x:x+handles.Param.pas2-1);
     imagesc(handles.sub1,'Parent',handles.im1)
     set(handles.im1,'Visible','off')
     
     x = handles.Posi(win2,1);
     y = handles.Posi(win2,2);
     handles.sub2=a(y:y+handles.Param.pas2-1,x:x+handles.Param.pas2-1);
     imagesc(handles.sub2,'Parent',handles.im2)
     set(handles.im2,'Visible','off')

     x = handles.Posi(win3,1);
     y = handles.Posi(win3,2);
     handles.sub3=a(y:y+handles.Param.pas2-1,x:x+handles.Param.pas2-1);
     imagesc(handles.sub3,'Parent',handles.im3)
     set(handles.im3,'Visible','off')
     handles.wins = [win1,win2,win3];
 
               
function handles = proporfft(handles)
     ic = handles.Param.pas2/2+1;
     jc = handles.Param.pas2/2+1;
     [I,J] = ndgrid(1:handles.Param.pas2,1:handles.Param.pas2);
     M = double((I-ic).^2+(J-jc).^2<=(handles.Param.cut)^2); 

    handles.fftsubcutsp1 = handles.fftsubcut1;
    handles.fftsubcutsp2 = handles.fftsubcut2;
    handles.fftsubcutsp3 = handles.fftsubcut3;
    
    
    reshaped = reshape(handles.fftsubcutsp1,[handles.Param.pas2*handles.Param.pas2 ,1]);         % Colonnes par colonnes donc l'indexage en 
    [~,index]= sort(reshaped ,'descend');
    nbel = handles.Param.propor * handles.Param.pas2*handles.Param.pas2;
    if nbel > handles.Param.pas2*handles.Param.pas2
        nbel =  handles.Param.pas2*handles.Param.pas2;
    else
    end
    
    handles.fftsubcutsp1(index(1:floor(nbel))) = 1;
    handles.fftsubcutsp1(index(floor(nbel)+1:end)) = 0;
    
 handles.fftsubcutsp1(M>0) = 1;

    spsubcutfilt = handles.fftsubcutsp1;
    B = false(handles.Param.pas2);
    B(spsubcutfilt>0)=1;
    se=strel('disk',handles.Param.strel);
    handles.fftsubcutspse1=imclose(B,se);
    spsubcutfilt =  handles.fftsubcutspse1;
    
    norm = sum(sum((spsubcutfilt).^2));
    [Yg,Xg]=meshgrid(1:(handles.Param.pas2),1:(handles.Param.pas2));
    xc = handles.Param.pas2/2;
    yc = handles.Param.pas2/2;
    Xo = Xg-xc; 
    Yo = Yg-yc;
    Y = (Yg-yc).^2; % Lignes
    X = (Xg-xc).^2; % Colonnes
    YY = sum(sum((spsubcutfilt).^2.*Y))/norm;
    XX = sum(sum((spsubcutfilt).^2.*X))/norm;
    XY = sum(sum((spsubcutfilt).^2.*Xo.*Yo))/norm;
    handles.M1 = [XX,XY;XY,YY];

     reshaped = reshape(handles.fftsubcutsp2,[handles.Param.pas2*handles.Param.pas2 ,1]);         % Colonnes par colonnes donc l'indexage en 
    [~,index]= sort(reshaped ,'descend');
    nbel = handles.Param.propor * handles.Param.pas2*handles.Param.pas2;
    if nbel > handles.Param.pas2*handles.Param.pas2
        nbel =  handles.Param.pas2*handles.Param.pas2;
    else
    end
    
    handles.fftsubcutsp2(index(1:floor(nbel))) = 1;
    handles.fftsubcutsp2(index(floor(nbel)+1:end)) = 0;
     handles.fftsubcutsp2(M>0) = 1;
 
    spsubcutfilt = handles.fftsubcutsp2;
    
    B = false(handles.Param.pas2);
    B(spsubcutfilt>0)=1;
    se=strel('disk',handles.Param.strel);
    handles.fftsubcutspse2=imclose(B,se);    
    spsubcutfilt =  handles.fftsubcutspse2;
    
    norm = sum(sum((spsubcutfilt).^2));
    [Yg,Xg]=meshgrid(1:(handles.Param.pas2),1:(handles.Param.pas2));
    xc = handles.Param.pas2/2;
    yc = handles.Param.pas2/2;
    Xo = Xg-xc; 
    Yo = Yg-yc;
    Y = (Yg-yc).^2; % Lignes
    X = (Xg-xc).^2; % Colonnes
    YY = sum(sum((spsubcutfilt).^2.*Y))/norm;
    XX = sum(sum((spsubcutfilt).^2.*X))/norm;
    XY = sum(sum((spsubcutfilt).^2.*Xo.*Yo))/norm;
    handles.M2 = [XX,XY;XY,YY];
    
    reshaped = reshape(handles.fftsubcutsp3,[handles.Param.pas2*handles.Param.pas2,1]);         % Colonnes par colonnes donc l'indexage en 
    [~,index]= sort(reshaped ,'descend');
    nbel = handles.Param.propor * handles.Param.pas2*handles.Param.pas2;
    if nbel > handles.Param.pas2*handles.Param.pas2
        nbel =  handles.Param.pas2*handles.Param.pas2;
    else
    end
    
    handles.fftsubcutsp3(index(1:floor(nbel))) = 1;
    handles.fftsubcutsp3(index(floor(nbel)+1:end)) = 0;
    handles.fftsubcutsp3(M>0) = 1;      

    spsubcutfilt = handles.fftsubcutsp3;
    B = false(handles.Param.pas2);
    B(spsubcutfilt>0)=1;
    se=strel('disk',handles.Param.strel);
    handles.fftsubcutspse3=imclose(B,se);  
    spsubcutfilt =  handles.fftsubcutspse3;
    
    norm = sum(sum((spsubcutfilt).^2));
    [Yg,Xg]=meshgrid(1:(handles.Param.pas2),1:(handles.Param.pas2));
    xc = handles.Param.pas2/2;
    yc = handles.Param.pas2/2;
    Xo = Xg-xc; 
    Yo = Yg-yc;
    Y = (Yg-yc).^2; % Lignes
    X = (Xg-xc).^2; % Colonnes
    YY = sum(sum((spsubcutfilt).^2.*Y))/norm;
    XX = sum(sum((spsubcutfilt).^2.*X))/norm;
    XY = sum(sum((spsubcutfilt).^2.*Xo.*Yo))/norm;
    handles.M3 = [XX,XY;XY,YY];
    
function handles = affichproporfft(handles)
        
   imagesc(handles.fftsubcutsp1,'Parent',handles.fft1propor)
   axes(handles.fft1propor);
   zoom(handles.Param.zoom)
 
   set(handles.fft1propor,'Visible','off')

   imagesc(handles.fftsubcutsp2,'Parent',handles.fft2propor)
    axes(handles.fft2propor);
   zoom(handles.Param.zoom)
   set(handles.fft2propor,'Visible','off')

   imagesc(handles.fftsubcutsp3,'Parent',handles.fft3propor)
       axes(handles.fft3propor);
   zoom(handles.Param.zoom)

   set(handles.fft3propor,'Visible','off')

   
function handles = trace_defo(handles)

    xo =  handles.Param.pas2/2;
    yo =  handles.Param.pas2/2;
        
    [~,E] = eig(handles.M1);
    L1 = 2*sqrt(E(2,2)/2);              
    L2 = 2*sqrt(E(1,1)/2);
    lgE = -1/2*(log(L1/L2));
    Dev = handles.M1-1/2*trace(handles.M1)*eye(2);
    [V,~] = eig(Dev);
    ang=atan(V(1,2)/V(2,2));
    handles.angS = ang;
    handles.S = lgE;

    amp = handles.S;
    ang = handles.angS;

    x1 = xo+cos(ang+pi/2)*amp*handles.Param.scale;
    y1 = yo+sin(ang+pi/2)*amp*handles.Param.scale;
    x2 = xo-cos(ang+pi/2)*amp*handles.Param.scale;
    y2 = yo-sin(ang+pi/2)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub1,'Parent',handles.im1)
    hold(handles.im1,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im1);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im1);
    set(handles.im1,'Visible','off')
    hold(handles.im1,'off');
    
    
    [~,E] = eig(handles.M2);
    L1 = 2*sqrt(E(2,2)/2);              
    L2 = 2*sqrt(E(1,1)/2);
    lgE = -1/2*(log(L1/L2));
    Dev = handles.M2-1/2*trace(handles.M2)*eye(2);
    [V,~] = eig(Dev);
    ang=atan(V(1,2)/V(2,2));
    angS = ang;
    S = lgE;

    amp = S;
    ang = angS;

    x1 = xo+cos(ang+pi/2)*amp*handles.Param.scale;
    y1 = yo+sin(ang+pi/2)*amp*handles.Param.scale;
    x2 = xo-cos(ang+pi/2)*amp*handles.Param.scale;
    y2 = yo-sin(ang+pi/2)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub2,'Parent',handles.im2)
    hold(handles.im2,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im2);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im2);
    set(handles.im2,'Visible','off')
    hold(handles.im2,'off');
    
     
    [~,E] = eig(handles.M3);
    L1 = 2*sqrt(E(2,2)/2);              
    L2 = 2*sqrt(E(1,1)/2);
    lgE = -1/2*(log(L1/L2));
    Dev = handles.M3-1/2*trace(handles.M3)*eye(2);
    [V,~] = eig(Dev);
    ang=atan(V(1,2)/V(2,2));
    angS = ang;
    S = lgE;

    amp = S;
    ang = angS;

    x1 = xo+cos(ang+pi/2)*amp*handles.Param.scale;
    y1 = yo+sin(ang+pi/2)*amp*handles.Param.scale;
    x2 = xo-cos(ang+pi/2)*amp*handles.Param.scale;
    y2 = yo-sin(ang+pi/2)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub3,'Parent',handles.im3)
    hold(handles.im3,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im3);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im3);
    set(handles.im3,'Visible','off')
    hold(handles.im3,'off');
    
function handles = ellipses_def(handles)
    abs_im_fft_w = handles.fftsubcut1;
    [yu,xu] = localMaximum_h(abs_im_fft_w,1,0,handles.Param.nbpoints);   
    
    [theta,rho]=cart2pol(xu-handles.Param.pas2/2,yu-handles.Param.pas2/2);
    handles.xu1=handles.Param.pas2./(2.*rho).*cos(theta);
    handles.yu1=handles.Param.pas2./(2.*rho).*sin(theta);
   
    set(handles.fft1cut,'Visible','on')
    imagesc(handles.fftsubcut1,'Parent',handles.fft1cut)
    hold(handles.fft1cut,'on');
    plot(handles.fft1cut,xu,yu,'Marker','+','Color','black','MarkerSize',10,'LineStyle','none')
   
    [~,~,ai,bi,phi,~]=ellipsefit(xu,yu);  
    hold(handles.fft1cut,'on');
    f_ellipseplotax(ai,bi,handles.Param.pas2/2+1,handles.Param.pas2/2+1,phi,'red',2,handles.fft1cut);
    handles.ai1 = ai;
    handles.bi1 =bi;
    handles.phi1 = phi;
    set(handles.fft1cut,'Visible','off')
    axes(handles.fft1cut);
    zoom(handles.Param.zoom)

    hold(handles.fft1cut,'off');
   clear xu yu
   
     abs_im_fft_w = handles.fftsubcut2;
    [yu,xu] = localMaximum_h(abs_im_fft_w,1,0,handles.Param.nbpoints);   
   
    set(handles.fft2cut,'Visible','on')
    imagesc(handles.fftsubcut2,'Parent',handles.fft2cut)
    hold(handles.fft2cut,'on');
    plot(handles.fft2cut,xu,yu,'Marker','+','Color','black','MarkerSize',10,'LineStyle','none')
   
    [~,~,ai,bi,phi,~]=ellipsefit(xu,yu);  
    hold(handles.fft2cut,'on');
    f_ellipseplotax(ai,bi,handles.Param.pas2/2+1,handles.Param.pas2/2+1,phi,'red',2,handles.fft2cut);
      handles.ai2 = ai;
    handles.bi2 =bi;
    handles.phi2 = phi;

    set(handles.fft2cut,'Visible','off')
   axes(handles.fft2cut);
    zoom(handles.Param.zoom)

    hold(handles.fft2cut,'off');
   clear xu yu
   
   
     abs_im_fft_w = handles.fftsubcut3;
    [yu,xu] = localMaximum_h(abs_im_fft_w,1,0,handles.Param.nbpoints);   
   
    set(handles.fft3cut,'Visible','on')
    imagesc(handles.fftsubcut3,'Parent',handles.fft3cut)
    hold(handles.fft3cut,'on');
    plot(handles.fft3cut,xu,yu,'Marker','+','Color','black','MarkerSize',10,'LineStyle','none')
   
    [~,~,ai,bi,phi,~]=ellipsefit(xu,yu);  
    hold(handles.fft3cut,'on');
    f_ellipseplotax(ai,bi,handles.Param.pas2/2+1,handles.Param.pas2/2+1,phi,'red',2,handles.fft3cut);
      handles.ai3 = ai;
    handles.bi3 =bi;
    handles.phi3 = phi;

    set(handles.fft3cut,'Visible','off')
   axes(handles.fft3cut);
    zoom(handles.Param.zoom);
    
    hold(handles.fft3cut,'off');

   
        
function handles = affichsize(handles)

    xo =  handles.Param.pas2/4;
    yo =  handles.Param.pas2/2;
        
    [~,E] = eig(handles.M1);
    L1 = sqrt(E(2,2)/2);              
    L2 = sqrt(E(1,1)/2);
    lgE = -1/2*(log(L1/L2));
    Dev = handles.M1-1/2*trace(handles.M1)*eye(2);
    [V,~] = eig(Dev);
    ang=atan(V(1,2)/V(2,2));
    handles.angS = ang;
    handles.S = lgE;

    amp = handles.S;
    ang = handles.angS;

    x1 = xo+cos(ang+pi/2)*amp*handles.Param.scale;
    y1 = yo+sin(ang+pi/2)*amp*handles.Param.scale;
    x2 = xo-cos(ang+pi/2)*amp*handles.Param.scale;
    y2 = yo-sin(ang+pi/2)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub1,'Parent',handles.im1)
    hold(handles.im1,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im1);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im1);
    hold(handles.im1,'on');
    f_ellipseplotax(handles.Param.pas2/2*1/handles.ai1,handles.Param.pas2/2*1/handles.bi1,xo+1/2*handles.Param.pas2,yo,handles.phi1,'red',2,handles.im1);
   
%     plot(handles.im1,handles.xu1+xo+1/2*handles.Param.pas2,handles.yu1+yo,'Marker','+','Color','yellow','MarkerSize',1,'LineStyle','none')     
%     Par =CircleFitByPratt([handles.xu1,handles.yu1]);  
%      hold(handles.im1,'on');
%      viscircles(handles.im1,Par(1,1:2),Par(1,3));
  %  f_ellipseplotax(ai,bi,xo+1/2*handles.Param.pas2,yo+1/4*handles.Param.pas2,phi,'yellow',handles.im1)
    
%     
    set(handles.im1,'Visible','off')
    hold(handles.im1,'off');
    
    
    [~,E] = eig(handles.M2);
    L1 = 2*sqrt(E(2,2)/2);              
    L2 = 2*sqrt(E(1,1)/2);
    lgE = -1/2*(log(L1/L2));
    Dev = handles.M2-1/2*trace(handles.M2)*eye(2);
    [V,~] = eig(Dev);
    ang=atan(V(1,2)/V(2,2));
    angS = ang;
    S = lgE;

    amp = S;
    ang = angS;

    x1 = xo+cos(ang+pi/2)*amp*handles.Param.scale;
    y1 = yo+sin(ang+pi/2)*amp*handles.Param.scale;
    x2 = xo-cos(ang+pi/2)*amp*handles.Param.scale;
    y2 = yo-sin(ang+pi/2)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    axes(handles.im2)
    imagesc(handles.sub2,'Parent',handles.im2)
    hold(handles.im2,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im2);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im2);
       hold(handles.im2,'on');
 
   f_ellipseplotax(handles.Param.pas2/2*1/handles.ai2,handles.Param.pas2/2*1/handles.bi2,xo+1/2*handles.Param.pas2,yo,handles.phi2,'red',2,handles.im2);
     hold(handles.im2,'off');
set(handles.im2,'Visible','off')
    
    [~,E] = eig(handles.M3);
    L1 = sqrt(E(2,2)/2);              
    L2 = sqrt(E(1,1)/2);
    lgE = -1/2*(log(L1/L2));
    Dev = handles.M3-1/2*trace(handles.M3)*eye(2);
    [V,~] = eig(Dev);
    ang=atan(V(1,2)/V(2,2));
    angS = ang;
    S = lgE;

    amp = S;
    ang = angS;

    x1 = xo+cos(ang+pi/2)*amp*handles.Param.scale;
    y1 = yo+sin(ang+pi/2)*amp*handles.Param.scale;
    x2 = xo-cos(ang+pi/2)*amp*handles.Param.scale;
    y2 = yo-sin(ang+pi/2)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub3,'Parent',handles.im3)
    hold(handles.im3,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im3);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im3);
    hold(handles.im3,'on');
    f_ellipseplotax(handles.Param.pas2/2*1/handles.ai3,handles.Param.pas2/2*1/handles.bi3,xo+1/2*handles.Param.pas2,yo,handles.phi3,'red',2,handles.im3);   
     set(handles.im3,'Visible','off')
    hold(handles.im3,'off');

    
    
function handles = affichsize_ellipsdef(handles)

    xo =  handles.Param.pas2/4;
    yo =  handles.Param.pas2/2;
        

    amp = 1/2*log(handles.bi1/handles.ai1);
    ang = handles.phi1+pi/2;

    x1 = xo+cos(ang)*amp*handles.Param.scale;
    y1 = yo+sin(ang)*amp*handles.Param.scale;
    x2 = xo-cos(ang)*amp*handles.Param.scale;
    y2 = yo-sin(ang)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub1,'Parent',handles.im1)
    hold(handles.im1,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im1);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im1);
    hold(handles.im1,'on');
    f_ellipseplotax(handles.Param.pas2/2*1/handles.ai1,handles.Param.pas2/2*1/handles.bi1,xo+1/2*handles.Param.pas2,yo,handles.phi1,'red',2,handles.im1);
       
    set(handles.im1,'Visible','off')
    hold(handles.im1,'off');
    
    
    amp = 1/2*log(handles.bi2/handles.ai2);
    ang = handles.phi2+pi/2;

    x1 = xo+cos(ang)*amp*handles.Param.scale;
    y1 = yo+sin(ang)*amp*handles.Param.scale;
    x2 = xo-cos(ang)*amp*handles.Param.scale;
    y2 = yo-sin(ang)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub2,'Parent',handles.im2)
    hold(handles.im2,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im2);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im2);
    hold(handles.im2,'on');
    f_ellipseplotax(handles.Param.pas2/2*1/handles.ai2,handles.Param.pas2/2*1/handles.bi2,xo+1/2*handles.Param.pas2,yo,handles.phi2,'red',2,handles.im2);
    set(handles.im2,'Visible','off')
    hold(handles.im2,'off');
    
        
    
    amp = 1/2*log(handles.bi3/handles.ai3);
    ang = handles.phi3+pi/2;

    x1 = xo+cos(ang)*amp*handles.Param.scale;
    y1 = yo+sin(ang)*amp*handles.Param.scale;
    x2 = xo-cos(ang)*amp*handles.Param.scale;
    y2 = yo-sin(ang)*amp*handles.Param.scale;

    X1 = [xo,x1];
    Y1 = [yo,y1];
    X2 = [xo,x2];
    Y2 = [yo,y2];
    
    imagesc(handles.sub3,'Parent',handles.im3)
    hold(handles.im3,'on');
    plot(X1,Y1,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im3);
    plot(X2,Y2,'Marker','.','LineStyle','-','Color','red','LineWidth',3,'Parent',handles.im3);
    hold(handles.im3,'on');
    f_ellipseplotax(handles.Param.pas2/2*1/handles.ai3,handles.Param.pas2/2*1/handles.bi3,xo+1/2*handles.Param.pas2,yo,handles.phi3,'red',2,handles.im3);  
    set(handles.im3,'Visible','off')
    hold(handles.im3,'off');
       axes(handles.fft3cut);
    zoom(handles.Param.zoom);
    
    hold(handles.fft3cut,'off');

    
    
function handles = filter(handles)
   imagesc(handles.fftsubcutspse1,'Parent',handles.filt1);
        axes(handles.filt1);
    zoom(handles.Param.zoom)
    set(handles.filt1,'Visible','off');
    
   imagesc(handles.fftsubcutspse2,'Parent',handles.filt2);
           axes(handles.filt2);
    zoom(handles.Param.zoom)

    set(handles.filt2,'Visible','off');
    
   imagesc(handles.fftsubcutspse3,'Parent',handles.filt3);
           axes(handles.filt3);
    zoom(handles.Param.zoom)
    set(handles.filt3,'Visible','off');

    
    


% --- Executes on button press in OKsize.
% function OKsize_Callback(hObject, eventdata, handles)
% % hObject    handle to OKsize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % hObject    handle to OKsize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% set(handles.OKsize,'String',handles.Param.nbpoints);
% guidata(hObject,handles);


% --- Executes on slider movement.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    handles.Param.zoom = get(hObject, 'Value');
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);

    handles = proporfft(handles);
    handles = affichproporfft(handles);
    handles = trace_defo(handles);
    
    set(handles.editzoom,'String',handles.Param.zoom);
    
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
    
    
    
        if get(handles.adv,'Value')
            handles = filter(handles);
        else
           
        end


guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
set(hObject,'Min',0,'Max',20,'SliderStep',[1/20,1/5],'Value',8);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function enterwindow_Callback(hObject, eventdata, handles)
% hObject    handle to enterwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enterwindow as text
%        str2double(get(hObject,'String')) returns contents of enterwindow as a double
    handles.Param.pas2 =  str2double(get(hObject,'String'));
    set(handles.window,'Value',handles.Param.pas2);
    if floor(handles.Param.pas2/2)==handles.Param.pas2/2
        handles.rec = handles.Param.pas2/2;
     else
         handles.rec = (handles.Param.pas2-1)/2;
    end

    handles = calculregion(handles);
    handles = affichim(handles);
    KB=handles.Param.pas2^2/16*9.53674e-7*8;  
    KBwhole = KB*handles.regl;
    set(handles.sizesub,'String',sprintf(['One image: ' num2str(floor(KBwhole)) ' MB']))
    handles = affichim(handles);
    handles = fftcalculus(handles);
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);
    handles = proporfft(handles);
    handles = affichproporfft(handles);
    handles = trace_defo(handles);

    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

    if handles.adv == 1
        handles = filter(handles);
    else
        set(get(handles.filt1,'children'),'visible','off')
        set(get(handles.filt2,'children'),'visible','off')
        set(get(handles.filt3,'children'),'visible','off')

    end




guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function enterwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enterwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editsigma_Callback(hObject, eventdata, handles)
% hObject    handle to editsigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editsigma as text
%        str2double(get(hObject,'String')) returns contents of editsigma as a double
    newval = str2double(get(hObject,'String'));
    handles.Param.sigma = newval;
    handles = fftcalculus(handles);
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);
    set(handles.sigma,'Value',newval);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

        
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function editsigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editsigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editcut_Callback(hObject, eventdata, handles)
% hObject    handle to editcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editcut as text
%        str2double(get(hObject,'String')) returns contents of editcut as a double
newval = str2double(get(hObject,'String'));
set(handles.cut,'Value',newval);
if  newval ~= handles.Param.cut
    set(hObject,'Value',newval);
    handles.Param.cut=newval;
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
      
else
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editpropor_Callback(hObject, eventdata, handles)
% hObject    handle to editpropor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editpropor as text
%        str2double(get(hObject,'String')) returns contents of editpropor as a double
    handles.Param.propor = str2double(get(hObject,'String'));
    set(handles.propor,'Value',handles.Param.propor*10^3);
   % handles.propor = handles.propor*10^-3;
    handles = proporfft(handles);
    handles = affichproporfft(handles);
    handles = trace_defo(handles);
    
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
    if get(handles.adv,'Value')
    handles = filter(handles);
    else
    end
   
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function editpropor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editpropor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editstrel_Callback(hObject, eventdata, handles)
% hObject    handle to editstrel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editstrel as text
%        str2double(get(hObject,'String')) returns contents of editstrel as a double
    newval = str2double(get(hObject,'String'));
    set(handles.strel,'Value',newval);
    if newval > 4 && newval ~= handles.Param.strel
        handles.Param.strel = newval;
        handles = proporfft(handles);
        handles = affichproporfft(handles);
        handles = trace_defo(handles);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
        handles = filter(handles);

    else
    end
guidata(hObject,handles);

% 


% --- Executes during object creation, after setting all properties.
function editstrel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editstrel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editsize_Callback(hObject, eventdata, handles)
% hObject    handle to editsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editsize as text
%        str2double(get(hObject,'String')) returns contents of editsize as a double

handles.Param.nbpoints = str2double(get(hObject,'String'));
set(handles.slidersize,'Value',handles.Param.nbpoints);
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end

% if handles.Param.sizeview==1 
%     handles = ellipses_def(handles);
%     handles = affichsize(handles);
% else
%     handles = affichcutfft(handles);
%     handles = trace_defo(handles);
% end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function editsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editzoom_Callback(hObject, eventdata, handles)
% hObject    handle to editzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editzoom as text
%        str2double(get(hObject,'String')) returns contents of editzoom as a double
    handles.Param.zoom = str2double(get(hObject, 'String'));
    handles = gaussfilt(handles);
    handles = cutfun(handles);
    handles = affichcutfft(handles);

    handles = proporfft(handles);
    handles = affichproporfft(handles);
    handles = trace_defo(handles);
    
    set(handles.zoom,'Value',handles.Param.zoom);
    
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')==0

        handles = ellipses_def(handles);
        handles = affichsize(handles);    
    else
    end
    if get(handles.sizes,'Value') && get(handles.anisotrpy_ellipse,'Value')
        handles = ellipses_def(handles);
        handles = affichsize_ellipsdef(handles);        
    end
        
        if get(handles.adv,'Value')
            handles = filter(handles);
        else
           
        end


guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function editzoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in anisotrpy_ellipse.
function anisotrpy_ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to anisotrpy_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anisotrpy_ellipse
handles.Param.ani_ellipse = get(hObject,'Value');

if get(hObject,'Value') &&  get(handles.sizes,'Value') 
    handles = ellipses_def(handles);
    handles = affichsize_ellipsdef(handles);
end
if get(hObject,'Value')==0 &&  get(handles.sizes,'Value') 
    handles = ellipses_def(handles);
    handles = affichsize(handles);
end
if get(hObject,'Value')==0 &&  get(handles.sizes,'Value')==0
    handles = affichcutfft(handles);
    handles = trace_defo(handles);
end

guidata(hObject,handles);


function nameani_ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to nameani_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nameani_ellipse as text
%        str2double(get(hObject,'String')) returns contents of nameani_ellipse as a double


% --- Executes during object creation, after setting all properties.
function nameani_ellipse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nameani_ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
