function varargout = Motor_Design_Rotor(varargin)
% MOTOR_DESIGN_ROTOR MATLAB code for Motor_Design_Rotor.fig
%      MOTOR_DESIGN_ROTOR, by itself, creates a new MOTOR_DESIGN_ROTOR or raises the existing
%      singleton*.
%
%      H = MOTOR_DESIGN_ROTOR returns the handle to a new MOTOR_DESIGN_ROTOR or the handle to
%      the existing singleton*.
%
%      MOTOR_DESIGN_ROTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTOR_DESIGN_ROTOR.M with the given input arguments.
%
%      MOTOR_DESIGN_ROTOR('Property','Value',...) creates a new MOTOR_DESIGN_ROTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Motor_Design_Rotor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Motor_Design_Rotor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Motor_Design_Rotor

% Last Modified by GUIDE v2.5 04-Sep-2017 01:33:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Motor_Design_Rotor_OpeningFcn, ...
                   'gui_OutputFcn',  @Motor_Design_Rotor_OutputFcn, ...
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


% --- Executes just before Motor_Design_Rotor is made visible.
function Motor_Design_Rotor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Motor_Design_Rotor (see VARARGIN)

% Choose default command line output for Motor_Design_Rotor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes1)
matlabImage = imread('Ranhura_Rotor.png');
image(matlabImage)
axis off
axis image

axes(handles.axes2)
matlabImage = imread('Ranhura_Estator_Rotor.png');
image(matlabImage)
axis off
axis image

axes(handles.axes3)
matlabImage = imread('Gaiola.png');
image(matlabImage)
axis off
axis image


% UIWAIT makes Motor_Design_Rotor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Motor_Design_Rotor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)


function edit1_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)


function edit2_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)


function edit3_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit4_Callback(hObject, eventdata, handles)


function edit4_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit5_Callback(hObject, eventdata, handles)


function edit5_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit6_Callback(hObject, eventdata, handles)


function edit6_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit7_Callback(hObject, eventdata, handles)


function edit7_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)


function edit8_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit9_Callback(hObject, eventdata, handles)


function edit9_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit10_Callback(hObject, eventdata, handles)


function edit10_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit11_Callback(hObject, eventdata, handles)


function edit11_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit12_Callback(hObject, eventdata, handles)


function edit12_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit13_Callback(hObject, eventdata, handles)


function edit13_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit14_Callback(hObject, eventdata, handles)


function edit14_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit15_Callback(hObject, eventdata, handles)


function edit15_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton2_Callback(hObject, eventdata, handles)

global raiz CP Sb Sa Nr ra g btte bttr hs2 hp2 bs2 hpeq rotor;
CP = str2num(get(handles.edit1,'String')); 
Sb = str2num(get(handles.edit2,'String'));
Nr = str2num(get(handles.edit3,'String'));
Sa = str2num(get(handles.edit4,'String'));
ra = str2num(get(handles.edit5,'String')); 
g = str2num(get(handles.edit7,'String'));
btte = str2num(get(handles.edit9,'String')); 
bttr = str2num(get(handles.edit18,'String'));
hs2 = str2num(get(handles.edit10,'String')); 
bs2 = str2num(get(handles.edit11,'String')); 
hp2 = str2num(get(handles.edit17,'String'));
hpeq = str2num(get(handles.edit15,'String')); 
rad1 = str2num(get(handles.edit20,'String')); 

rotor.slots = Nr;
rotor.gap = g;
rotor.rad2 = (rotor.rad2 - g) * 1e3;
rotor.rad1 = rad1 * 1e3;
rotor.magnet_poles = 2;

save_as = [get(handles.edit19,'String') '.mat'];
global re bs HC P Nep Ns Nea RIE HRa HRp hs1a phi_p phi_a pho_cond bsm hp1 bp1 hs1 NC;
global Kw1 alpha Bs Br wm lm l y urec Bpk_dente Bpk_coroa Bg  t h1 h2 D stator
cd(raiz);
cd programa;
Calculo_Geometrico;
cd(raiz);
cd modelos
save(save_as,'Ntotp','Ntota','LCabp','LCaba','Neep','Neea','LCtotp',...
    'LCtota','Ltotp','Ltota','A_fiop','A_fioa','Vtotp','Vtota','Peso_p','Peso_a',...
    'bs','HC','Nep','Nea','Ns','P','RIE','HRa','HRp','hs1a','phi_p','phi_a','pho_cond','bsm','hp1','hs1','bp1','NC',...
    'CP','Sb','Sa','Nr','ra','g','btte','bttr','hp2','bs2','hs2','hpeq', 'ref', 'Xxm', 'pho_rotor', 'pho_cond', 'MIA', ...
    'MIP', 'alpha', 'Bs', 'Br', 'wm', 'lm', 'l', 'y', 'urec', 're', 'rotor', 'stator');
cd(raiz);
cd interface;
close;
menu
movegui(menu,'center')

% global Ntotp Ntota LCabp LCaba Neep Neea LCtotp LCtota Ltotp Ltota A_fiop A_fioa Vtotp Vtota Peso_p Peso_a

function pushbutton3_Callback(hObject, eventdata, handles)
close
menu
movegui(menu,'center')


function edit17_Callback(hObject, eventdata, handles)


function edit17_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit18_Callback(hObject, eventdata, handles)


function edit18_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit19_Callback(hObject, eventdata, handles)


function edit19_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
