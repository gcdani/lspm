function varargout = LSPM_Plot_Manager(varargin)
% LSPM_PLOT_MANAGER MATLAB code for LSPM_Plot_Manager.fig
%      LSPM_PLOT_MANAGER, by itself, creates a new LSPM_PLOT_MANAGER or raises the existing
%      singleton*.
%
%      H = LSPM_PLOT_MANAGER returns the handle to a new LSPM_PLOT_MANAGER or the handle to
%      the existing singleton*.
%
%      LSPM_PLOT_MANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSPM_PLOT_MANAGER.M with the given input arguments.
%
%      LSPM_PLOT_MANAGER('Property','Value',...) creates a new LSPM_PLOT_MANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSPM_Plot_Manager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSPM_Plot_Manager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSPM_Plot_Manager

% Last Modified by GUIDE v2.5 12-Aug-2016 01:24:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSPM_Plot_Manager_OpeningFcn, ...
                   'gui_OutputFcn',  @LSPM_Plot_Manager_OutputFcn, ...
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


% --- Executes just before LSPM_Plot_Manager is made visible.
function LSPM_Plot_Manager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSPM_Plot_Manager (see VARARGIN)

% Choose default command line output for LSPM_Plot_Manager
handles.output = hObject;
set(hObject,'toolbar','figure');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LSPM_Plot_Manager wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSPM_Plot_Manager_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton2.
axes(handles.axes1);
global delta P_entrada T_frente Eff P_cu_estator_tras P_cu_estator_frente P_entreferro T_rel Wcur;

switch get(handles.popupmenu1,'Value') 
    case 1
        plot_option = 1;
    case 2
        plot_option = 2;
    case 3
        plot_option = 3;
    case 4
        plot_option = 4;
    case 5
        plot_option = 5;
    case 6
        plot_option = 6;
    case 7
        plot_option = 7;
    case 8
        plot_option = 8;
    case 9
        plot_option = 9;
    otherwise
end

if plot_option == 1
    
    plot(delta,P_entrada);
    grid on
    title('Potencia de entrada (Pot. elec no SPEED)');
    xlabel('Delta (graus)');
    ylabel('Potencia (W)');

elseif plot_option == 2

    plot(delta, T_frente);
    grid on
    title('Torque pra frente');
    xlabel('Delta (graus)');
    ylabel('Potencia(W)')
    
elseif plot_option == 3
    
    plot(delta,  Eff);
    grid on
    title('Torque pra tras');
    xlabel('Delta (graus)');
    ylabel('Potencia(Nm)')
    
elseif plot_option == 4
    
    plot(delta,P_cu_estator_frente);
    grid on
    title('Potencia Cu estator para frente');
    xlabel('Delta (graus)');
    ylabel('Potencia (W)');
    
elseif plot_option == 5
    
    plot(delta, P_cu_estator_tras);
    grid on
    title('Potencia Cu estator para tras');
    xlabel('Delta (graus)');
    ylabel('Potencia (W)');
    
elseif plot_option == 6
    
    plot(delta,Wcur);
    grid on
    title('Potencia no condutor do rotor');
    xlabel('Delta (graus)');
    ylabel('Potencia (W)');
    
elseif plot_option == 7
    
    plot(delta,P_entreferro);
    grid on
    title('Pot. entreferro');
    xlabel('Delta (graus)');
    ylabel('Potencia (W)');
    
elseif plot_option == 8
    
    plot(delta,T_rel);
    grid on
    title('Eficiencia');
    xlabel('Delta (graus)');
    ylabel('Eficiencia');
    
elseif plot_option == 9
    
    plot( Eff,T_frente);
    grid on
    title('Eff');
    xlabel('T-frente');
    ylabel('Eff');
    
end

guidata(hObject, handles); 

function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset')
guidata(hObject, handles); %updates the handles


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all
global raiz
cd(raiz)