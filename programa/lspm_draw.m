%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a LSPM motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw motor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stator
% -----
% stator ring
mi_drawarc(0, stat_outer_rad, 0, -stat_outer_rad, 180, 1);
mi_drawarc(0, -stat_outer_rad, 0, stat_outer_rad, 180, 1);
mi_selectnode(0, stat_outer_rad);
mi_selectnode(0, -stat_outer_rad);
mi_selectarcsegment(stat_outer_rad, 0);
mi_selectarcsegment(-stat_outer_rad, 0);
mi_setgroup(group_stator);
% stator teeth
mi_drawpolyline(stat_slot_shape_r);
mi_drawpolyline(stat_slot_shape_l);
for i = 1:length(stat_slot_shape_r)
    mi_selectnode(stat_slot_shape_r(i,:));
    mi_selectnode(stat_slot_shape_l(i,:));
end
mi_setgroup(group_stator);
for i = 1:length(stat_slot_shape_r)
    mi_selectnode(stat_slot_shape_r(i,:));
    mi_selectnode(stat_slot_shape_l(i,:));
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = 1:length(stat_slot_shape_r) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
end
mi_setgroup(group_stator);
for i = 1:length(stat_slot_shape_r) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 1);
mi_addarc(stat_base_n_x, stat_base_n_y, stat_slot_shape_r(1,1), stat_slot_shape_r(1,2), stat_slot_base_angle * 2, 1);
mi_selectarcsegment((stat_base_n_x + stat_slot_shape_r(1,1)) / 2, (stat_base_n_y + stat_slot_shape_r(1,2)) / 2);
mi_setgroup(group_stator);
mi_selectarcsegment((stat_base_n_x + stat_slot_shape_r(1,1)) / 2, (stat_base_n_y + stat_slot_shape_r(1,2)) / 2);
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 3);

if SSlotBot == 0
    mult_arc_rot = 12;
elseif SSlotBot == 1
    mult_arc_rot = 0.1;
end

mi_addarc(stat_slot_shape_r(end,1), stat_slot_shape_r(end,2), stat_slot_shape_l(end,1), stat_slot_shape_l(end,2), stat_slot_base_angle * mult_arc_rot, 1);
mi_selectarcsegment((stat_slot_shape_l(end,1) + stat_slot_shape_r(end,1)) / 2, (stat_slot_shape_l(end,2) + stat_slot_shape_r(end,2)) / 2);
mi_setgroup(group_stator);
mi_selectarcsegment((stat_slot_shape_l(end,1) + stat_slot_shape_r(end,1)) / 2, (stat_slot_shape_l(end,2) + stat_slot_shape_r(end,2)) / 2);
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 3);
    
mi_clearselected();

% Windings
mi_zoomnatural();
mi_drawpolyline(stat_main_wind_cont_r);
mi_drawpolyline(stat_main_wind_cont_l);
for i = 1:length(stat_main_wind_cont_r)
    mi_selectnode(stat_main_wind_cont_r(i,:));
    mi_selectnode(stat_main_wind_cont_l(i,:));
end
mi_setgroup(group_stator);
for i = 1:length(stat_main_wind_cont_r)
    mi_selectnode(stat_main_wind_cont_r(i,:));
    mi_selectnode(stat_main_wind_cont_l(i,:));
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = 1:length(stat_main_wind_cont_r) - 1
    mi_selectsegment((stat_main_wind_cont_r(i,:) + stat_main_wind_cont_r(i+1,:)) / 2);
    mi_selectsegment((stat_main_wind_cont_l(i,:) + stat_main_wind_cont_l(i+1,:)) / 2);
end
mi_setgroup(group_stator);
for i = 1:length(stat_main_wind_cont_r) - 1
    mi_selectsegment((stat_main_wind_cont_r(i,:) + stat_main_wind_cont_r(i+1,:)) / 2);
    mi_selectsegment((stat_main_wind_cont_l(i,:) + stat_main_wind_cont_l(i+1,:)) / 2);
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 1);

mi_drawpolyline(stat_aux_wind_cont_r);
mi_drawpolyline(stat_aux_wind_cont_l);
for i = 1:length(stat_aux_wind_cont_r)
    mi_selectnode(stat_aux_wind_cont_r(i,:));
    mi_selectnode(stat_aux_wind_cont_l(i,:));
end
mi_setgroup(group_stator);
for i = 1:length(stat_aux_wind_cont_r)
    mi_selectnode(stat_aux_wind_cont_r(i,:));
    mi_selectnode(stat_aux_wind_cont_l(i,:));
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = 1:length(stat_aux_wind_cont_r) - 1
    mi_selectsegment((stat_aux_wind_cont_r(i,:) + stat_aux_wind_cont_r(i+1,:)) / 2);
    mi_selectsegment((stat_aux_wind_cont_l(i,:) + stat_aux_wind_cont_l(i+1,:)) / 2);
end
mi_setgroup(group_stator);
for i = 1:length(stat_aux_wind_cont_r) - 1
    mi_selectsegment((stat_aux_wind_cont_r(i,:) + stat_aux_wind_cont_r(i+1,:)) / 2);
    mi_selectsegment((stat_aux_wind_cont_l(i,:) + stat_aux_wind_cont_l(i+1,:)) / 2);
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 1);


% Rotor
% ------
% inner circle
mi_drawarc(0, rot_out_rad, 0, -rot_out_rad, 180, 1);
mi_drawarc(0, -rot_out_rad, 0, rot_out_rad, 180, 1);
mi_drawarc(0, rot_inner_rad, 0, -rot_inner_rad, 180, 1);
mi_drawarc(0, -rot_inner_rad, 0, rot_inner_rad, 180, 1);
mi_selectnode(0, rot_out_rad);
mi_selectnode(0, -rot_out_rad);
mi_selectnode(0, rot_inner_rad);
mi_selectnode(0, -rot_inner_rad);
mi_selectarcsegment(rot_out_rad, 0);
mi_selectarcsegment(-rot_out_rad, 0);
mi_selectarcsegment(rot_inner_rad, 0);
mi_selectarcsegment(-rot_inner_rad, 0);
mi_setgroup(group_rotor);
% slots
mi_drawpolyline(rot_slot_shape_r);
mi_drawpolyline(rot_slot_shape_l);
for i = 1:length(rot_slot_shape_r)
    mi_selectnode(rot_slot_shape_r(i,:));
    mi_selectnode(rot_slot_shape_l(i,:));
end
mi_setgroup(group_rotor);
for i = 1:length(rot_slot_shape_r)
    mi_selectnode(rot_slot_shape_r(i,:));
    mi_selectnode(rot_slot_shape_l(i,:));
end
mi_copyrotate2(0, 0, 360 / rot_nof_slots, rot_nof_slots - 1, 0);
for i = 1:length(rot_slot_shape_r) - 1
    mi_selectsegment((rot_slot_shape_r(i,:) + rot_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((rot_slot_shape_l(i,:) + rot_slot_shape_l(i+1,:)) / 2);
end
mi_setgroup(group_rotor);
for i = 1:length(rot_slot_shape_r) - 1
    mi_selectsegment((rot_slot_shape_r(i,:) + rot_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((rot_slot_shape_l(i,:) + rot_slot_shape_l(i+1,:)) / 2);
end
mi_copyrotate2(0, 0, 360 / rot_nof_slots, rot_nof_slots - 1, 1);
mi_addarc(rot_slot_shape_r(1,1), rot_slot_shape_r(1,2), rot_base_n_x, rot_base_n_y, rot_slot_base_angle * 12, 1);
mi_selectarcsegment((rot_base_n_x + rot_slot_shape_r(1,1)) / 2, (rot_base_n_y + rot_slot_shape_r(1,2)) / 2);
mi_setgroup(group_rotor);
mi_selectarcsegment((rot_base_n_x + rot_slot_shape_r(1,1)) / 2, (rot_base_n_y + rot_slot_shape_r(1,2)) / 2);
mi_copyrotate2(0, 0, 360 / rot_nof_slots, rot_nof_slots - 1, 3);

% magnets slots
mi_drawpolygon(rot_mag_shape_s_r);
mi_drawpolygon(rot_mag_shape_s_l);
mi_drawpolygon(rot_mag_shape_i_r);
mi_drawpolygon(rot_mag_shape_i_l);
for i = 1:length(rot_mag_shape_s_r)
    mi_selectnode(rot_mag_shape_s_r(i,:));
    mi_selectnode(rot_mag_shape_s_l(i,:));
    mi_selectnode(rot_mag_shape_i_r(i,:));
    mi_selectnode(rot_mag_shape_i_l(i,:));
end
mi_setgroup(group_rotor);

for i = 1:length(rot_mag_shape_s_r) - 1
    mi_selectsegment((rot_mag_shape_s_r(i,:) + rot_mag_shape_s_r(i+1,:)) / 2);
    mi_selectsegment((rot_mag_shape_s_l(i,:) + rot_mag_shape_s_l(i+1,:)) / 2);
    mi_selectsegment((rot_mag_shape_i_r(i,:) + rot_mag_shape_i_r(i+1,:)) / 2);
    mi_selectsegment((rot_mag_shape_i_l(i,:) + rot_mag_shape_i_l(i+1,:)) / 2);
end
mi_setgroup(group_rotor);

% magnets
mi_drawpolygon(rot_mag_s_r);
mi_drawpolygon(rot_mag_s_l);
mi_drawpolygon(rot_mag_i_r);
mi_drawpolygon(rot_mag_i_l);
for i = 1:length(rot_mag_s_r)
    mi_selectnode(rot_mag_s_r(i,:));
    mi_selectnode(rot_mag_s_l(i,:));
    mi_selectnode(rot_mag_i_r(i,:));
    mi_selectnode(rot_mag_i_l(i,:));
end
mi_setgroup(group_rotor);

for i = 1:length(rot_mag_s_r) - 1
    mi_selectsegment((rot_mag_s_r(i,:) + rot_mag_s_r(i+1,:)) / 2);
    mi_selectsegment((rot_mag_s_l(i,:) + rot_mag_s_l(i+1,:)) / 2);
    mi_selectsegment((rot_mag_i_r(i,:) + rot_mag_i_r(i+1,:)) / 2);
    mi_selectsegment((rot_mag_i_l(i,:) + rot_mag_i_l(i+1,:)) / 2);
end
mi_setgroup(group_rotor);

% Boundary around motor
% ---------------------
mi_zoomnatural();
mi_drawarc(0, stat_outer_rad * 3, 0, -stat_outer_rad * 3, 180, 1);
mi_drawarc(0, -stat_outer_rad * 3, 0, stat_outer_rad * 3, 180, 1);

% Block labels
% ------------
% Air
mi_addblocklabel(0, 0);
mi_selectlabel(0, 0);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_addblocklabel(0, 2 * stat_outer_rad);
mi_selectlabel(0, 2 * stat_outer_rad);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_addblocklabel(0, rot_out_rad + gap / 2);
mi_selectlabel(0,  rot_out_rad + gap / 2);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_clearselected();
% rotor iron
mi_addblocklabel(0, (rot_inner_rad + rot_base_rad) / 2);
mi_selectlabel(0, (rot_inner_rad + rot_out_rad) / 2);
mi_setgroup(group_rotor);
mi_selectlabel(0, (rot_inner_rad + rot_out_rad) / 2);
mi_setblockprop(mat_rot, 1, 0, '<none>', 0, group_rotor, 0)
mi_clearselected();
% stator iron
mi_addblocklabel(0, (stat_outer_rad + stat_inner_rad) / 2);
mi_selectlabel(0, (stat_outer_rad + stat_inner_rad) / 2);
mi_setgroup(group_stator);
mi_selectlabel(0, (stat_outer_rad + stat_inner_rad) / 2);
mi_setblockprop(mat_stat, 1, 0, '<none>', 0, group_stator, 0)
mi_clearselected();
% Magnets
mi_addblocklabel([mag_shape_center_up_x  mag_shape_center_up_y]);
mi_addblocklabel([mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_addblocklabel([mag_shape_center_up_x  -mag_shape_center_up_y]);
mi_addblocklabel([mag_shape_center_bot_x -mag_shape_center_bot_y]);
mi_addblocklabel([-mag_shape_center_up_x  mag_shape_center_up_y]);
mi_addblocklabel([-mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_addblocklabel(-[mag_shape_center_up_x  mag_shape_center_up_y]);
mi_addblocklabel(-[mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_selectlabel([mag_shape_center_up_x  mag_shape_center_up_y]);
mi_selectlabel([mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_selectlabel([mag_shape_center_up_x  -mag_shape_center_up_y]);
mi_selectlabel([mag_shape_center_bot_x -mag_shape_center_bot_y]);
mi_selectlabel([-mag_shape_center_up_x  mag_shape_center_up_y]);
mi_selectlabel([-mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_selectlabel(-[mag_shape_center_up_x  mag_shape_center_up_y]);
mi_selectlabel(-[mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0);
mi_selectlabel([mag_shape_center_up_x  mag_shape_center_up_y]);
mi_selectlabel([mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_selectlabel([mag_shape_center_up_x  -mag_shape_center_up_y]);
mi_selectlabel([mag_shape_center_bot_x -mag_shape_center_bot_y]);
mi_selectlabel([-mag_shape_center_up_x  mag_shape_center_up_y]);
mi_selectlabel([-mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_selectlabel(-[mag_shape_center_up_x  mag_shape_center_up_y]);
mi_selectlabel(-[mag_shape_center_bot_x mag_shape_center_bot_y]);
mi_setgroup(group_rotor);


mi_addblocklabel([mag_center_x  mag_center]);
mi_addblocklabel([mag_center_x  -mag_center]);
mi_addblocklabel([-mag_center_x  mag_center]);
mi_addblocklabel(-[mag_center_x  mag_center]);
mi_selectlabel([mag_center_x  mag_center]);
mi_selectlabel([mag_center_x  -mag_center]);
mi_selectlabel([-mag_center_x  mag_center]);
mi_selectlabel(-[mag_center_x  mag_center]);
mi_setgroup(group_rotor);

dir = 90 - gama;

mi_selectlabel([mag_center_x  mag_center]);
mi_setblockprop(mat_magnet, 1, 0, '<none>', dir, group_rotor, 0);
mi_clearselected();

mi_selectlabel([mag_center_x  -mag_center]);
mi_setblockprop(mat_magnet, 1, 0, '<none>', -dir, group_rotor, 0);
mi_clearselected();

mi_selectlabel([-mag_center_x  mag_center]);
mi_setblockprop(mat_magnet, 1, 0, '<none>', 180 - dir, group_rotor, 0);
mi_clearselected();

mi_selectlabel(-[mag_center_x  mag_center]);
mi_setblockprop(mat_magnet, 1, 0, '<none>', 180 + dir, group_rotor, 0);
mi_clearselected();

% Windings
for i = 0:stat_nof_slots-1
    main_wind_angle_p = 2 * pi / stat_nof_slots * i + stat_main_wind_angle;
    main_wind_angle_n = 2 * pi / stat_nof_slots * i - stat_main_wind_angle;
    aux_wind_angle_p = 2 * pi / stat_nof_slots * i - stat_aux_wind_angle;
    aux_wind_angle_n = 2 * pi / stat_nof_slots * i + stat_aux_wind_angle;
    mi_addblocklabel([ stat_main_wind_rad * cos(main_wind_angle_p)
                       stat_main_wind_rad * sin(main_wind_angle_p)]);
    mi_addblocklabel([ stat_main_wind_rad * cos(main_wind_angle_n)
                       stat_main_wind_rad * sin(main_wind_angle_n)]);
    mi_addblocklabel([ stat_aux_wind_rad * cos(aux_wind_angle_p)
                       stat_aux_wind_rad * sin(aux_wind_angle_p)]);
    mi_addblocklabel([ stat_aux_wind_rad * cos(aux_wind_angle_n)
                       stat_aux_wind_rad * sin(aux_wind_angle_n)]);
    mi_selectlabel([ stat_main_wind_rad * cos(main_wind_angle_p)
                       stat_main_wind_rad * sin(main_wind_angle_p)]);
    mi_selectlabel([ stat_main_wind_rad * cos(main_wind_angle_n)
                       stat_main_wind_rad * sin(main_wind_angle_n)]);
    mi_selectlabel([ stat_aux_wind_rad * cos(aux_wind_angle_p)
                       stat_aux_wind_rad * sin(aux_wind_angle_p)]);
    mi_selectlabel([ stat_aux_wind_rad * cos(aux_wind_angle_n)
                       stat_aux_wind_rad * sin(aux_wind_angle_n)]);               
    mi_setgroup(group_stator);
    mi_selectlabel([ stat_main_wind_rad * cos(main_wind_angle_p)
                       stat_main_wind_rad * sin(main_wind_angle_p)]);
    mi_setblockprop(mat_wind, 1, 0, 'L1+', 0, group_stator, stat_nof_wdg)
    mi_clearselected();
    mi_selectlabel([ stat_main_wind_rad * cos(main_wind_angle_n)
                       stat_main_wind_rad * sin(main_wind_angle_n)]);
    mi_setblockprop(mat_wind, 1, 0, 'L1-', 0, group_stator, stat_nof_wdg)
    mi_clearselected();
    
    mi_selectlabel([ stat_aux_wind_rad * cos(aux_wind_angle_p)
                       stat_aux_wind_rad * sin(aux_wind_angle_p)]);
    mi_setblockprop(mat_wind, 1, 0, 'L2+', 0, group_stator, stat_nof_wdg)
    mi_clearselected();
    mi_selectlabel([ stat_aux_wind_rad * cos(aux_wind_angle_n)
                       stat_aux_wind_rad * sin(aux_wind_angle_n)]);
    mi_setblockprop(mat_wind, 1, 0, 'L2-', 0, group_stator, stat_nof_wdg)
    mi_clearselected();

end
mi_clearselected();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom view to current problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi_zoomnatural();
