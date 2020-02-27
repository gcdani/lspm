function varargout = LSPM_Plot_Manager_Assincrono(varargin)
% LSPM_PLOT_MANAGER_ASSINCRONO MATLAB code for LSPM_Plot_Manager_Assincrono.fig
%      LSPM_PLOT_MANAGER_ASSINCRONO, by itself, creates a new LSPM_PLOT_MANAGER_ASSINCRONO or raises the existing
%      singleton*.
%
%      H = LSPM_PLOT_MANAGER_ASSINCRONO returns the handle to a new LSPM_PLOT_MANAGER_ASSINCRONO or the handle to
%      the existing singleton*.
%
%      LSPM_PLOT_MANAGER_ASSINCRONO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSPM_PLOT_MANAGER_ASSINCRONO.M with the given input arguments.
%
%      LSPM_PLOT_MANAGER_ASSINCRONO('Property','Value',...) creates a new LSPM_PLOT_MANAGER_ASSINCRONO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSPM_Plot_Manager_Assincrono_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSPM_Plot_Manager_Assincrono_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSPM_Plot_Manager_Assincrono

% Last Modified by GUIDE v2.5 25-Aug-2016 13:03:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSPM_Plot_Manager_Assincrono_OpeningFcn, ...
                   'gui_OutputFcn',  @LSPM_Plot_Manager_Assincrono_OutputFcn, ...
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


% --- Executes just before LSPM_Plot_Manager_Assincrono is made visible.
function LSPM_Plot_Manager_Assincrono_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSPM_Plot_Manager_Assincrono (see VARARGIN)

% Choose default command line output for LSPM_Plot_Manager_Assincrono
handles.output = hObject;
set(hObject,'toolbar','figure');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LSPM_Plot_Manager_Assincrono wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSPM_Plot_Manager_Assincrono_OutputFcn(hObject, eventdata, handles) 
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
warning off
global Rot_sinc Tind Tmag Tf Tt Tqdi Tqqi Vid Viq Tmotor

switch get(handles.popupmenu1,'Value') 
    case 1
        plot(0:Rot_sinc,Tmotor,'k');
        hold on
        plot(0:Rot_sinc,Tind,'--b');
        plot(0:Rot_sinc,Tmag,'--r');
        legend('Tmotor','Tind','Tmag');
        grid on
        title('Motor Torque');
        xlabel('Speed (RPM)');
        ylabel('Torque (kgf.cm)');
        axis([-inf,inf,0,inf])
        hold off
    case 2
        plot(0:Rot_sinc,abs(Tf));
        grid on
        title('Direct Torque');
        xlabel('Speed (RPM)');
        ylabel('Torque (kgf.cm)');
        axis([-inf,inf,0,inf])
    case 3
        plot(0:Rot_sinc,abs(Tt));
        grid on
        title('Inverse Torque');
        xlabel('Speed (RPM)');
        ylabel('Torque (kgf.cm)');
        axis([-inf,inf,0,inf])
    case 4
        plot(0:Rot_sinc,Tind);
        grid on
        title('Motor Torque - no Magnet Torque');
        xlabel('Speed (RPM)');
        ylabel('Torque (kgf.cm)');
        axis([-inf,inf,0,inf])
    case 5
        plot(0:Rot_sinc,Tmag,'k');
        hold on
        plot(0:Rot_sinc,abs(Tqdi),'--r');
        plot(0:Rot_sinc,abs(Tqqi),'--g');
        legend('Total Resistive Torque','Resistive Torque Direct','Resistive Torque Quadrature');
        grid on
        title('Magnet Torque');
        xlabel('Speed (RPM)');
        ylabel('Torque (kgf.cm)');
        axis([-inf,inf,0,inf])
        hold off
    case 6
        plot(0:Rot_sinc,Vid,'r');
        hold on
        plot(0:Rot_sinc,Viq,'--b');
        grid on
        legend('FEM - Direct Axis','FEM - Quadrature Axis');
        title('Induced Voltage');
        xlabel('Speed (RPM)');
        ylabel('Induced Voltage (V)');
        axis([-inf,inf,0,inf])
        hold off
    otherwise
end

warning on

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


function pushbutton4_Callback(hObject, eventdata, handles)
close
LSPM_Sincrono


function pushbutton5_Callback(hObject, eventdata, handles)
[name, path] = uiputfile('dados_plot_aasinc.txt');
path_file = fullfile(path, name);
fid = fopen(path_file,'wt');
global Rot_sinc Tind Tmag Tf Tt Tqdi Tqqi Vid Viq Tmotor;
switch get(handles.popupmenu1,'Value') 
    case 1
        fprintf(fid, 'RPM            T_motor        T_ind          T_mag\n');
        fprintf(fid,'%f\t%f\t%f\t%f\n',[0:Rot_sinc;Tmotor; Tind; Tmag]);
    case 2
        fprintf(fid, 'RPM            T_frente\n');
        fprintf(fid,'%f\t%f\n',[0:Rot_sinc;abs(Tf)]);
    case 3
        fprintf(fid, 'RPM            T_tras\n');
        fprintf(fid,'%f\t%f\n',[0:Rot_sinc;abs(Tt)]);
    case 4
        fprintf(fid, 'RPM            T_ind\n');
        fprintf(fid,'%f\t%f\n',[0:Rot_sinc;Tind]);
    case 5
        fprintf(fid, 'RPM            T_ag           T_qdi          T_qqi\n');
        fprintf(fid,'%f\t%f\t%f\t%f\n',[0:Rot_sinc;Tmag; abs(Tqdi); abs(Tqqi)]);
    case 6
        fprintf(fid, 'RPM            Vid            Viq\n');
        fprintf(fid,'%f\t%f\t%f\n',[0:Rot_sinc;Vid; Viq]);
    otherwise
end

fclose(fid);
