function varargout = LSPM_Plot_Manager_Sincrono(varargin)
% LSPM_PLOT_MANAGER_SINCRONO MATLAB code for LSPM_Plot_Manager_Sincrono.fig
%      LSPM_PLOT_MANAGER_SINCRONO, by itself, creates a new LSPM_PLOT_MANAGER_SINCRONO or raises the existing
%      singleton*.
%
%      H = LSPM_PLOT_MANAGER_SINCRONO returns the handle to a new LSPM_PLOT_MANAGER_SINCRONO or the handle to
%      the existing singleton*.
%
%      LSPM_PLOT_MANAGER_SINCRONO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSPM_PLOT_MANAGER_SINCRONO.M with the given input arguments.
%
%      LSPM_PLOT_MANAGER_SINCRONO('Property','Value',...) creates a new LSPM_PLOT_MANAGER_SINCRONO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSPM_Plot_Manager_Sincrono_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSPM_Plot_Manager_Sincrono_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSPM_Plot_Manager_Sincrono

% Last Modified by GUIDE v2.5 04-Sep-2017 01:34:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSPM_Plot_Manager_Sincrono_OpeningFcn, ...
                   'gui_OutputFcn',  @LSPM_Plot_Manager_Sincrono_OutputFcn, ...
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


% --- Executes just before LSPM_Plot_Manager_Sincrono is made visible.
function LSPM_Plot_Manager_Sincrono_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSPM_Plot_Manager_Sincrono (see VARARGIN)

% Choose default command line output for LSPM_Plot_Manager_Sincrono
handles.output = hObject;
set(hObject,'toolbar','figure');
% Update handles structure
guidata(hObject, handles);
cla(handles.axes1,'reset')
% UIWAIT makes LSPM_Plot_Manager_Sincrono wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSPM_Plot_Manager_Sincrono_OutputFcn(hObject, eventdata, handles) 
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
global WCu_main WCu_aux PF delta Pelec Tps Eff WCu1 WCu2 Pshaft T_rel WCuR T_total Tns T_pma;


switch get(handles.popupmenu1,'Value') 
    case 1
        plot(delta,Pelec);
        grid on
        title('Potencia de Entrada');
        xlabel('Delta (Graus)');
        ylabel('Potencia (W)');
        hold off
    case 2
        plot(delta, Tps,'k');
        grid on
        title('Torque para Frente');
        xlabel('Delta (Graus)');
        ylabel('Torque (Kgf.cm)');
        hold on
        plot(delta, T_rel.*10.2,'--b');
        plot(delta, T_pma.*10.2,'--r');
        legend('Torque para frente', 'Torque relutancia', 'Torque de alinhamento');
        hold off
    case 3
        plot(delta, Tns.*10.2);
        grid on
        title('Torque para Tras');
        xlabel('Delta (Graus)');
        ylabel('Torque (Kgf.cm)');
    case 4
        plot(delta,WCu1);
        grid on
        title('Potencia Cu Estator para Frente');
        xlabel('Delta (Graus)');
        ylabel('Potencia (W)');
    case 5
        plot(delta, WCu2);
        grid on
        title('Potencia Cu Estator para Tras');
        xlabel('Delta (Graus)');
        ylabel('Potencia (W)');
        axis([-inf,inf,0,inf])
    case 6
        plot(delta,WCuR);
        grid on
        title('Potencia no Condutor do Rotor');
        xlabel('Delta (Graus)');
        ylabel('Potencia (W)');
    case 7
        plot(delta,Pshaft);
        grid on
        title('Potencia de Saida');
        xlabel('Delta (Graus)');
        ylabel('Potencia (W)');
    case 8
        plot(delta,Eff);
        grid on
        title('Efficiencia x Delta');
        xlabel('Delta (Graus)');
        ylabel('Efficiencia (%)');
    case 9
        plot(T_total, Eff);
        grid on
        title('Efficiencia x Torque de Saida');
        xlabel('Torque (kgf.cm)');
        ylabel('Efficiencia (%)');
    case 10
        plot(delta, T_rel.*10.2);
        grid on
        title('Torque de Relutancia');
        xlabel('Delta (Graus)');
        ylabel('Torque (Kgf.cm)');
    case 11
        plot(delta, T_pma.*10.2);
        grid on
        title('Torque de Alinhamento');
        xlabel('Delta (Graus)');
        ylabel('Torque (Kgf.cm)');
    case 12
        plot(delta, T_total,'k');
        grid on
        title('Torque Total');
        xlabel('Delta (Graus)');
        ylabel('Torque (Kgf.cm)');
        hold on
        plot(delta, Tps.*10.2,'--b');
        plot(delta, Tns.*10.2,'--r');
        legend('Torque total', 'Torque para frente', 'Torque para tras');
        hold off
    case 13
        plot(delta, PF);
        grid on
        title('Fator de Potencia');
        xlabel('Delta (Graus)');
        ylabel('Fator de Potencia');
    case 14
        plot(Eff, WCu_main, '--b');
        grid on
        title('Perdas no Cobre no Enrolamento Principal e Auxiliar pela Eficiencia');
        xlabel('Eficiencia (%)');
        ylabel('Perdas no Cobre (W)');
        hold on
        plot(Eff, WCu_aux, '--r');
        legend('Perda no Enrolamento Principal', 'Perdas no Enrolamento Auxiliar');
        hold off
    otherwise
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


function pushbutton4_Callback(hObject, eventdata, handles)
[name, path] = uiputfile('dados_plot_sinc.txt');
path_file = fullfile(path, name);
fid = fopen(path_file,'wt');
global delta_carga P_in T_frente Eff P_cu_tras P_cu_frente P_out T_tras Wcur T_total T_rel T_al;
switch get(handles.popupmenu1,'Value') 
    case 1
        fprintf(fid, 'Delta          Potencia de Entrada \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;P_in]);
    case 2
        fprintf(fid, 'Delta          Torque Frente  T_rel          T_al\n');
        fprintf(fid,'%f\t%f\t%f\t%f\n',[delta_carga;T_frente; T_rel; T_al]);
    case 3
        fprintf(fid, 'Delta          Torque para tras \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;T_tras]);
    case 4
        fprintf(fid, 'Delta          P_cu_estator_frente \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;P_cu_frente]);
    case 5
        fprintf(fid, 'Delta          P_cu_estator_tras \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;P_cu_tras]);
    case 6
        fprintf(fid, 'Delta          Wcur \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;Wcur]);
    case 7
        fprintf(fid, 'Delta          Potencia de Saida \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;P_out]);
    case 8
        fprintf(fid, 'Delta          Efficiencia \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;Eff]);
    case 9
        fprintf(fid, 'T_total        Efficiencia \n');
        fprintf(fid,'%f\t%f\n',[T_total;Eff]);
    case 10
        fprintf(fid, 'Delta          T_rel \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;T_rel]);
    case 11
        fprintf(fid, 'Delta          T_al \n');
        fprintf(fid,'%f\t%f\n',[delta_carga;T_al]);
    case 12
        fprintf(fid, 'Delta          T_total        T_frente       T_tras\n');
        fprintf(fid,'%f\t%f\t%f\t%f\n',[delta_carga;T_total;T_frente;T_tras]);
    otherwise
end

fclose(fid);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
