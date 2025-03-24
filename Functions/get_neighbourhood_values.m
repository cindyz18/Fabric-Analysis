function [neighbourhood, neighbourhood_location]= get_neighbourhood_values(P, contact_list)
% This function extracts the 26 neighboring voxel values around each contact point in a 3D matrix.
% Inputs:
%   P: 3D matrix representing labeled particles, where each voxel is assigned a particle label.
%   contact_list: Nx4 matrix where each row contains a contact label and its (X, Y, Z) coordinates.
%
% Outputs:
%   neighbourhood: Nx26 matrix where each row contains the values of the 26 neighboring voxels for each contact point.
%   neighbourhood_location: Nx3 matrix storing the (X, Y, Z) coordinates of each contact point.

[a1, a2, a3] = size(P);
counter = 1;
% Preallocate for maximum possible size
neighbourhood_location = zeros(size(contact_list,1), 3);
neighbourhood = zeros(size(contact_list,1), 26);

for idx = 1:size(contact_list,1)
    i = contact_list(idx, 2);
    j = contact_list(idx, 3);
    k = contact_list(idx, 4);
    % Ensure we are not at the boundary to avoid out-of-bounds error
    if i > 1 && i < a1 && j > 1 && j < a2 && k > 1 && k < a3
        % Extract the 3x3x3 neighborhood - getting the 27 values - 26 neighbours and R_voxel value in a matrix
        neighbors = P(i-1:i+1, j-1:j+1, k-1:k+1);

        % Reshape and remove the central voxel value
        % Reshaping the matrix into a single row array
        contact_neighbors = reshape(neighbors, 1, []);
        contact_neighbors(14) = [];

        % Store the location and neighbors
        neighbourhood_location(counter, :) = [i, j, k];
        neighbourhood(counter, :) = contact_neighbors;

        counter = counter + 1;
    end
end

% Remove unused preallocated rows
neighbourhood_location = neighbourhood_location(1:counter-1, :);
neighbourhood = neighbourhood(1:counter-1, :);

end