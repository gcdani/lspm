%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a LSPM motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants for defining motor dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stator
stat_nof_poles = stator.poles; % 2
mot_thickness = stator.thick; % 10
stat_nof_wdg = stator.main_wdg; % 50
stat_nof_slots = stator.slots; % 24
SSlotBot = stator.SlotBot; % 0 = round; 1 = square
stat_outer_rad = stator.rad; % 42.44
w1s = stator.w1s; % 0.5
h1s = stator.h1s; % 0.25
w2s = stator.w2s; % 2
h2s = stator.h2s; % 0.5
w3s = stator.w3s; % 4
h3s = stator.h3s; % 6
filSB = 0; % nao feito
filSO = 0; % nao feito

% Rotor
rot_nof_mag_per_pole = rotor.magnet_poles; % 2
rot_nof_slots = rotor.slots; % 28
gap = rotor.gap; % 0.35
RMI = rotor.rmi; % 17.1122
Inset = rotor.inset; % 1
MagWid = rotor.magwid; % 14
MalotWid = rotor.malotwid; % 21.3
LM = rotor.lm; % 2.5
CWeb = rotor.cweb; % 0.5
Vtrap = rotor.vtrap; % 117 degrees
rot_out_rad = rotor.rad2; % 29.52
rot_inner_rad = rotor.rad1; % 11
rot_base_rad = rotor.base; % 20
rot_head_angle = rotor.head_angle; % 45
SD_R = rotor.sd; % 6
R_Bridge = rotor.bridge; % 0.25
TW_R = rotor.tw; % 2.68
rot_slot_width = 2 * pi * rot_out_rad / rot_nof_slots - TW_R;
stat_inner_rad = rot_out_rad + gap;

% Magnet


% rotor
rot_slot_shape_r = [
                    TW_R / 2    rot_out_rad - R_Bridge - SD_R;
                    TW_R / 2    rot_out_rad - R_Bridge - 0.5 * rot_slot_width * sin(rot_head_angle * pi / 180);
                    (TW_R + rot_slot_width) / 2     rot_out_rad - R_Bridge
                    ];

rot_slot_shape_l = [
                    -TW_R / 2    rot_out_rad - R_Bridge - SD_R;
                    -TW_R / 2    rot_out_rad - R_Bridge - 0.5 * rot_slot_width * sin(rot_head_angle * pi / 180);
                    -(TW_R + rot_slot_width) / 2     rot_out_rad - R_Bridge
                    ];

% stator
stat_slot_shape_r = [ w1s / 2     stat_inner_rad;
w1s / 2     stat_inner_rad + h1s;
w2s / 2     stat_inner_rad + h1s + h2s;
w3s / 2     stat_inner_rad + h1s + h2s + h3s];

stat_slot_shape_l = [
-w1s / 2     stat_inner_rad;
-w1s / 2     stat_inner_rad + h1s;
-w2s / 2     stat_inner_rad + h1s + h2s;
-w3s / 2     stat_inner_rad + h1s + h2s + h3s];                
                
mat_air     = 'Air';
mat_stat    = 'Vanadium Permedur';
mat_rot     = 'Vanadium Permedur';
mat_magnet  = 'NdFeB 52 MGOe';
%mat_magnet  = 'NdFeB 40 MGOe';
mat_wind    = '1mm';

current = 5;

% Groups for stator and rotor
group_stator = 1;
group_rotor = 2;
group_air = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate needed variables for motor dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotor edge
rot_head_width = rot_slot_width + TW_R;
rot_head_edge = sqrt((rot_out_rad)^2 - (rot_head_width/2)^2);
% rotor base edge
rot_base_edge = rot_base_rad * cos(asin(rot_slot_width / 2 / rot_base_rad));
% rotor base angle
rot_slot_edge_angle = atan(rot_slot_shape_r(1,1) / rot_slot_shape_r(1,2)) *180 / pi;
rot_slot_base_angle = 360 / rot_nof_slots / 2 - rot_slot_edge_angle;
% Stator winding contour
stat_aux_wind_cont_r = [ 0                  stat_inner_rad + h1s + h2s;
                        w2s / 2             stat_inner_rad + h1s + h2s;
                        (w2s + w3s) / 4     stat_inner_rad + h1s + h2s + h3s / 2;
                        0                   stat_inner_rad + h1s + h2s + h3s / 2];
                    
stat_aux_wind_cont_l = [  -w2s / 2          stat_inner_rad + h1s + h2s;
                        0                   stat_inner_rad + h1s + h2s;
                        0                   stat_inner_rad + h1s + h2s + h3s / 2;
                        -(w2s + w3s) / 4    stat_inner_rad + h1s + h2s + h3s / 2];                    
                    
stat_main_wind_cont_r = [0    stat_inner_rad + h1s + h2s + h3s / 2;
                        (w2s + w3s) / 4     stat_inner_rad + h1s + h2s + h3s / 2;
                        w3s / 2             stat_inner_rad + h1s + h2s + h3s;
                        0                   stat_inner_rad + h1s + h2s + h3s];  
                    
stat_main_wind_cont_l = [ -(w2s + w3s) / 4   stat_inner_rad + h1s + h2s + h3s / 2;
                        0                    stat_inner_rad + h1s + h2s + h3s / 2;
                        0                    stat_inner_rad + h1s + h2s + h3s;
                        -w3s / 2             stat_inner_rad + h1s + h2s + h3s];  

% Stator winding radius
stat_main_wind_rad = sqrt(stat_main_wind_cont_r(3,1) ^ 2 + ( stat_inner_rad + h1s + h2s + h3s * 0.75 ) ^ 2);                   
stat_aux_wind_rad =  sqrt(stat_aux_wind_cont_r(3,1) ^ 2 + ( stat_inner_rad + h1s + h2s + h3s / 4 ) ^ 2);                   
% Stator winding angle
stat_main_wind_angle = atan( 0.5 * stat_main_wind_cont_r(3,1) / stat_main_wind_cont_r(3,2) );
stat_aux_wind_angle = atan( 0.5 * stat_aux_wind_cont_r(3,1) / stat_aux_wind_cont_r(3,2) );
% magnet shape
gama = Vtrap * 0.5 * pi / 180;
mag_x = RMI / sin(gama) - MalotWid * cos(gama);
mag_y = MalotWid * sin(gama);
mag_end = (RMI + LM) / sin(gama);
rot_mag_shape_s_l = [ 
                 RMI / sin(gama) - CWeb * 0.5 / tan(gama)     CWeb * 0.5;
                 mag_end                                      CWeb * 0.5;
                 mag_end + (0.5 * CWeb - mag_y) / tan(gama)   mag_y;
                 mag_x                                        mag_y;
                 RMI / sin(gama) - CWeb * 0.5 / tan(gama)     CWeb * 0.5];
             
rot_mag_shape_s_r = rot_mag_shape_s_l;             
rot_mag_shape_s_r(:,1) = -rot_mag_shape_s_r(:,1);

rot_mag_shape_i_l = rot_mag_shape_s_l;             
rot_mag_shape_i_l(:,2) = -rot_mag_shape_i_l(:,2);

rot_mag_shape_i_r = - rot_mag_shape_s_l;

mag_center = RMI / tan(gama) + Inset;
mag_center_2 = (RMI + LM) / tan(gama) + Inset;
mag_bot = mag_center - MagWid / 2;
mag_bot_2 = mag_center_2 - MagWid / 2;
mag_up = mag_center + MagWid / 2;
mag_up_2 = mag_center_2 + MagWid / 2;
mag_bot_x = RMI / sin(gama) - mag_bot * cos(gama);
mag_bot_x_2 = (RMI + LM) / sin(gama) - mag_bot_2 * cos(gama);
mag_bot_y = mag_bot * sin(gama);
mag_bot_y_2 = mag_bot_2 * sin(gama);
mag_up_x = RMI / sin(gama) - mag_up * cos(gama);
mag_up_x_2 = (RMI + LM) / sin(gama) - mag_up_2 * cos(gama);
mag_up_y = mag_up * sin(gama);
mag_up_y_2 = mag_up_2 * sin(gama);

rot_mag_s_l = [ 
                 mag_bot_x     mag_bot_y;
                 mag_bot_x_2   mag_bot_y_2;
                 mag_up_x_2    mag_up_y_2;
                 mag_up_x      mag_up_y;
                 mag_bot_x     mag_bot_y];
             
rot_mag_s_r = rot_mag_s_l;             
rot_mag_s_r(:,1) = -rot_mag_s_r(:,1);

rot_mag_i_l = rot_mag_s_l;             
rot_mag_i_l(:,2) = -rot_mag_i_l(:,2);

rot_mag_i_r = - rot_mag_s_l;

mag_shape_center_up_y = ( mag_y + mag_up_y ) / 2;
mag_shape_center_up_x = ( mag_x + mag_end + (0.5 * CWeb - mag_y) / tan(gama) ) / 2;

mag_shape_center_bot_y = ( CWeb / 2 + mag_bot_y ) / 2;
mag_shape_center_bot_x = ( RMI / sin(gama) - CWeb * 0.5 / tan(gama) + mag_end ) / 2;

mag_center_x = RMI / sin(gama) - mag_center / tan(gama);
mag_center_x_2 = (RMI + LM) / sin(gama) - mag_center_2 / tan(gama);
mag_center_x = ( mag_center_x + mag_center_x_2 ) / 2;

% position of next slot
rot_base_n_x = rot_base_rad * sin((rot_slot_edge_angle + 2 * rot_slot_base_angle) / 180 * pi);
rot_base_n_y = rot_base_rad * cos((rot_slot_edge_angle + 2 * rot_slot_base_angle) / 180 * pi);

stat_slot_edge_angle = atan(stat_slot_shape_r(1,1) / stat_slot_shape_r(1,2)) *180 / pi;
stat_slot_base_angle = 360 / stat_nof_slots / 2 - stat_slot_edge_angle;
stat_base_n_x = stat_inner_rad * sin((stat_slot_edge_angle + 2 * stat_slot_base_angle) / 180 * pi);
stat_base_n_y = stat_inner_rad * cos((stat_slot_edge_angle + 2 * stat_slot_base_angle) / 180 * pi);