function varargout = LSPM_Sincrono(varargin)
%LSPM_SINCRONO M-file for LSPM_Sincrono.fig
%      LSPM_SINCRONO, by itself, creates a new LSPM_SINCRONO or raises the existing
%      singleton*.
%
%      H = LSPM_SINCRONO returns the handle to a new LSPM_SINCRONO or the handle to
%      the existing singleton*.
%
%      LSPM_SINCRONO('Property','Value',...) creates a new LSPM_SINCRONO using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to LSPM_Sincrono_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LSPM_SINCRONO('CALLBACK') and LSPM_SINCRONO('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LSPM_SINCRONO.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSPM_Sincrono

% Last Modified by GUIDE v2.5 04-Sep-2017 01:36:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSPM_Sincrono_OpeningFcn, ...
                   'gui_OutputFcn',  @LSPM_Sincrono_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before LSPM_Sincrono is made visible.
function LSPM_Sincrono_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for LSPM_Sincrono
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global raiz

cd ..
cd modelos
aux1 = dir('*.mat');
for i=1:1:size(aux1)
    casos_exist{i}=aux1(i).name;
end
if exist('casos_exist') ~= 0 %se n?o encontrou nenhum .mat desabilita o carregar, sen?o preencha a lista
      set(handles.listbox1,'String',casos_exist);
end
cd ..
cd interface


% UIWAIT makes LSPM_Sincrono wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSPM_Sincrono_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function freqRede_Callback(hObject, eventdata, handles)
% hObject    handle to freqRede (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqRede as text
%        str2double(get(hObject,'String')) returns contents of freqRede as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
set(hObject,'String','0') 
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function freqRede_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqRede (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Vm_rms_Callback(hObject, eventdata, handles)
% hObject    handle to Vm_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vm_rms as text
%        str2double(get(hObject,'String')) returns contents of Vm_rms as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
set(hObject,'String','0') 
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Vm_rms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vm_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Cstart_Callback(hObject, eventdata, handles)
% hObject    handle to Cstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Cstart as text
%        str2double(get(hObject,'String')) returns contents of Cstart as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
set(hObject,'String','0') 
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Cstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI Bsat w raiz T f Cap Vlin uo Crun;
f = str2num(get(handles.freqRede,'String')); 
Vlin = str2num(get(handles.Vm_rms,'String')); 
Cap = str2num(get(handles.Cstart,'String'));
Crun = str2num(get(handles.Crun,'String'));
T = str2num(get(handles.edit5,'String'));
select   = get(handles.listbox1,'Value'); % N?mero da consulta selecionada
caso_sel = get(handles.listbox1,'String');% Nomes das consultas
dados_nome = strtok(caso_sel{select},'.');
cd(raiz)
cd modelos
load(caso_sel{select})% Carrega vari?veis 
cd(raiz)
cd programa

w = 2*pi*f;
CI = 1;
Bsat = 1.1;
uo = 4*pi*1e-7;
Calculo_Resistencias;
Calculo_Indutancias;
Calculo_Inducao;

Parte_Sincrona2


function edit5_Callback(hObject, eventdata, handles)


function edit5_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton2_Callback(hObject, eventdata, handles)
warning off
global w raiz T f Cap Vlin CI Bsat uo Crun
f = str2num(get(handles.freqRede,'String')); 
Vlin = str2num(get(handles.Vm_rms,'String')); 
Cap = str2num(get(handles.Cstart,'String'));
Crun = str2num(get(handles.Crun,'String'));
T = str2num(get(handles.edit5,'String'));
select   = get(handles.listbox1,'Value'); % N?mero da consulta selecionada
caso_sel = get(handles.listbox1,'String');% Nomes das consultas
dados_nome = strtok(caso_sel{select},'.');
cd(raiz)
cd modelos
load(caso_sel{select})% Carrega vari?veis 
cd(raiz)
cd programa

w = 2*pi*f;
CI = 1;
Bsat = 1.1;
uo = 4*pi*1e-7;
Calculo_Resistencias;
Calculo_Indutancias;
Calculo_Inducao;

Parte_Assincrona
warning on

function pushbutton3_Callback(hObject, eventdata, handles)

close all
menu;
movegui(LSPM_Sincrono,'center')


function Crun_Callback(hObject, eventdata, handles)


function Crun_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
