%% Porosity Mapper Script
clear;clc;
addpath('./../Functions','./../Data')
%% Define Inputs
% Define the input file name
file_name = 'UG_0pc_sand_binary_v1';
rev = 103;

% Define output file location
output_file_name = 'UG_0pc_sand_binary_porosity';
output_file_location = '../Results/';
file_name_prefix = strcat(output_file_location, output_file_name);

%% Calculating REV value
tic; % Start timer
disp('Loading variables')
% Load the .mat file containing the 3D matrix
disp(['Loading file: ', file_name]);
load(strcat(file_name, '.mat'));

% Extract matrix using dynamic naming (avoid changing "Avizo_" and "_mat")
eval(strcat('M0=Avizo_',file_name,'_mat;'));
M0 = squeeze(M0); % Remove singleton dimensions

disp("Mapping porosity")
[a1, a2, a3] = size(M0);
ROI = rev*rev*rev;
stored = calculate_porosity_values(M0, rev);
porosity=stored(:,1)./ROI;

% Compute the average porosity
porosity_avg = sum(porosity) / ((a1-rev) * (a2-rev) * (a3-rev));
disp(['The average porosity is ', num2str(porosity_avg)])

time_taken = toc;
disp(['Time taken by movsum approach: ', num2str(time_taken), ' seconds']);

tic; % Start timer
B = map_porosity_to_matrix(porosity, a1, a2, a3, rev);
time_taken = toc;
disp(['Time taken to map porosity: ', num2str(time_taken), ' seconds']);


%% Save Results
% save(strcat(file_name_prefix,'.mat'),'porosity','-v7.3')
% disp("Porosity matrix saved")

% Save matrix
save(strcat(file_name_prefix,'_B.mat'),'B','-v7.3')
disp("B matrix saved")

