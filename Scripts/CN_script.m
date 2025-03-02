% Coordination Number Script
clc;close all;clear all
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
eval(strcat('labeled=Avizo_',sand_file_name,'_mat;'))
labeled=squeeze(labeled);

contacts = contacts(201:400,201:400,201:400);
labeled = labeled(201:400,201:400,201:400); 

%% Coordination Number
[particle_list, contact_list, labeled_modified] = get_contact_list(labeled, contacts);
[neighbourhood, neighbourhood_location]= get_neighbourhood_values(labeled_modified, contact_list);
[value_counts] = get_coordination_number(neighbourhood, labeled);
% disp("Coordination number : [Voxel value, Repetitions]")
% value_counts
[centroids, part_coords_cell] = get_particle_centroids(labeled_modified, particle_list);




%% Plotting
% Plotting for large data sizes will take time to run, this is for
% visualisation purposes only, do not run section if plotting is not
% required
figure('Name',file_name);

% Plotting the contacts
plot3(contact_list(:,2),contact_list(:,3),contact_list(:,4),'.','MarkerSize',20,'Color','yellow');       
hold on

% Plotting the particles
for c = 1:size(part_coords_cell,1)
    if ~isempty(part_coords_cell{c})
        part_coords = part_coords_cell{c};
        k = boundary(part_coords(:, 1), part_coords(:, 2), part_coords(:, 3));
        trisurf(k, part_coords(:, 1), part_coords(:, 2), part_coords(:, 3), ...
            'Facecolor', 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none');
    end
end
colorbar;
hold on

% Plotting the centroids of particles as points
plot3(centroids(:, 3), centroids(:, 4), centroids(:, 5), '.', 'MarkerSize', 20, 'Color', 'green');

saveas(gcf, strcat(file_name, 'script', '.png'));