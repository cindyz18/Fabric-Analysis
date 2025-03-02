function [value_counts] = get_coordination_number(neighbours, P)
% get_coordination_number computes the coordination number of particle contacts
% Inputs:
%   neighbours: NxM matrix where each row represents neighboring particle labels
%   P: 3D matrix defining the labelled particles
% Outputs:
%   value_counts: Nx2 matrix containing unique voxel values and their counts


% Preallocate the output matrix M3 with zeros
[a1, ~] = size(neighbours);
M3 = zeros(a1, 2);
ConVoxVal = max(P ,[], "all" ) + 1;
disp(['ConVoxVal: ', num2str(ConVoxVal)])

% Replace zeros with ConVoxVal and get unique values row-wise
neighbours(neighbours == 0) = ConVoxVal;

for i = 1:a1
    unique_neighbours = unique(neighbours(i, :)); % Get unique values in the row
    ss1 = numel(unique_neighbours);                % Number of unique values

    if ss1 == 3
        M3(i, :) = unique_neighbours(1:2);         % Assign first two unique values
    end
end

% Remove rows with all zeros
unique_rows = unique(M3(M3(:, 1) ~= 0, :), 'rows');

disp("Unique sets of combinations")
[voxel_value,~,ic] = unique(unique_rows);
a_counts_opt = accumarray(ic,1);
value_counts = [voxel_value, a_counts_opt];

end