%% REV FINDER Optimisation Script
% REV FINDER Optimisation Script
% Author: Cindy Zhu
% Date: 23/02/205=025
% Version: 1.0
%
% Description:
% This script calculates porosity values for various Representative 
% Elementary Volume (REV) sizes and positions within a 3D matrix loaded 
% from a .mat file. The results are saved to an Excel file and visualized 
% in a plot.
%
% Inputs:
% - A .mat file containing a 3D matrix (e.g., 'UG_0pc_sand_binary_v1.mat').
% - REV sizes (default: 3:10:250).
% - Positions for porosity calculations (4 predefined coordinates).
%
% Outputs:
% - An Excel file ('REVfinder.xlsx') with porosity data.
% - A PNG plot visualizing porosity vs. REV size.
%
% Usage:
% Simply run the script in MATLAB. Ensure the .mat file is in the same 
% directory or update the file path accordingly.
%
% Dependencies:
% - Requires 'calculate_porosity_values.m' function.
% - MATLAB toolbox for 'xlswrite' and 'movsum' functions.
%
% Notes:
% - Avoid modifying the "Avizo_" and "_mat" parts of the filename.
% - Ensure the Excel file is not open during writing to prevent errors.

clear;clc;
addpath('./../Functions','./../Data')
%% Define Inputs
% Define the input file name
file_name = 'UG_0pc_sand_binary_v1';
% Define REV sizes and positions
REV = 3:10:250; % Array of REV sizes, this will need to be chnaged depending on Matrix size
pos = [350 350 505; 350 650 505; 650 350 505; 650 650 505]; % Coordinate positions in matrix to conduct REV

%% Calculating REV value
tic; % Start timer
disp('Loading variables')
% Load the .mat file containing the 3D matrix
disp(['Loading file: ', file_name]);
load(strcat(file_name, '.mat'));

% Extract matrix using dynamic naming (avoid changing "Avizo_" and "_mat")
eval(strcat('M0=Avizo_',file_name,'_mat;'));
M0 = squeeze(M0); % Remove singleton dimensions

% Define REV sizes and positions
[y1, y2] = size(REV);
[x1, x2] = size(pos);

% Preallocate matrix to store porosity values
M2 = zeros(x1 * y2, 5);
counter = 1;

disp('Calculating REV values for each slice of the matrix')
% Loop through each REV size
for REV_index = 1:y2
    len = REV(REV_index) - 1;
    for ii = 1:x1
        % Slice the matrix based on position and REV size
        M1 = M0(...
            pos(ii,1)-(len/2):pos(ii,1)+(len/2),...
            pos(ii,2)-(len/2):pos(ii,2)+(len/2),...
            pos(ii,3)-(len/2):pos(ii,3)+(len/2)...
        );

        % Calculate porosity using the helper function
        [a1, a2, a3] = size(M1);
        L = 1; % Porosity region of interest size
        stored = calculate_porosity_values(M1, L);

        % Compute the average porosity
        porosity_avg = sum(stored) / ((a1-L) * (a2-L) * (a3-L));

        % Store REV, position, and porosity
        M2(counter, :) = [REV(REV_index), pos(ii,1), pos(ii,2), pos(ii,3), porosity_avg];
        counter = counter + 1;
    end
end

%% Saving results
disp('Writing results to file')
% Write results to Excel
col_header = {'REV size', 'Loc_x', 'Loc_y', 'Loc_z', 'Porosity'};
xlswrite('REVfinder.xlsx', M2, file_name, 'B2');
xlswrite('REVfinder.xlsx', col_header, file_name, 'B1');

% Display execution time
time_taken = toc;
disp(['Time taken: ', num2str(time_taken), ' seconds']);

%% Plotting Results
disp('Plotting results')
% Plotting porosity vs REV Size for each unique coordinate
unique_rows = unique(M2(:, 2:4), 'rows');
figure('Name',file_name);
hold on;
for row_idx = 1:size(unique_rows, 1)
    row = unique_rows(row_idx, :);
    indices = ismember(M2(:, 2:4), row, 'rows');
    plot(M2(indices, 1), M2(indices, 5), '-o', 'DisplayName', sprintf('(%d, %d, %d)', row));
end

% Format the plot
xlabel('REV Size');
ylabel('Porosity');
title('Porosity vs REV Size for Different Coordinates');
legend('Location', 'best');
grid on;
hold off;

% Save the plot as a PNG file
saveas(gcf, strcat(file_name, '.png'));
