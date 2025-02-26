function [value_counts] = get_coordination_number(neighbours, P)
% get_coordination_number computes the coordination number of particle contacts
% Inputs:
%   neighbours: 3D matrix defining the labelled particles
% Outputs:
%   value_counts: Matrix with unique voxel values and their counts

%%% Getting the unique values of the rows and storing them

% Start timer for performance tracking
tic;

% Preallocate the output matrix M3 with zeros
[a1, ~] = size(neighbours);
M3 = zeros(a1, 2);
ConVoxVal = max(P ,[], "all" ) + 1;
disp(['ConVoxVal: ', num2str(ConVoxVal)])

% Replace zeros with ConVoxVal and get unique values row-wise
M2_modified_opt = neighbours;
M2_modified_opt(M2_modified_opt == 0) = ConVoxVal;

for i = 1:a1
    M_temp = unique(M2_modified_opt(i, :)); % Get unique values in the row
    ss1 = numel(M_temp);                % Number of unique values

    if ss1 == 3
        M3(i, :) = M_temp(1:2);         % Assign first two unique values
    end
end

% Remove rows with all zeros
M3_uniq_opt = unique(M3(M3(:, 1) ~= 0, :), 'rows');

disp("Unique sets of combinations")
[voxel_value,~,ic] = unique(M3_uniq_opt);
a_counts_opt = accumarray(ic,1);
value_counts = [voxel_value, a_counts_opt];

end