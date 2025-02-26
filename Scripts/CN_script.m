% Coordination Number Script
% clc;close all;clear all
addpath('./../Functions','./../Data')

%% Define Inputs
% Define the input file name
file_name = 'UG_0pc_sand_binary_v1';
% 'UG_0pc_sand_contact_binary_v1' 'UG_0pc_sand_separated_labelled_v1'
contact_file_name = 'UG_0pc_sand_contact_binary_v1';
sand_file_name = 'UG_0pc_sand_separated_labelled_v1';

load(strcat(contact_file_name, '.mat')); %Evaluate the first 3 lines - then copy the generated matrix name into line 4.
eval(strcat('contacts=Avizo_',contact_file_name,'_mat;'))  %Do not remove or change the first 6 and last 4 characters - "Avizo_" and "_mat"
contacts=squeeze(contacts);

load(strcat(sand_file_name, '.mat'));
eval(strcat('labelled=Avizo_',sand_file_name,'_mat;'))
labelled=squeeze(labelled);

contacts = contacts(201:400,201:400,201:400);
labelled = labelled(201:400,201:400,201:400); 

%% Coordination Number
[particle_list, contact_list, labelled_modified] = get_contact_list(labelled, contacts);
[neighbourhood, neighbourhood_location]= get_neighbourhood_values(labelled_modified, contact_list);
[value_counts] = get_coordination_number(neighbourhood,labelled);
% disp("Coordination number : [Voxel value, Repetitions]")
% value_counts

%% Plotting
figure('Name',file_name);
plot3(contact_list(:,2),contact_list(:,3),contact_list(:,4),'.','MarkerSize',20,'Color','yellow');       
hold on