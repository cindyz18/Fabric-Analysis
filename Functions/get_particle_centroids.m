function [centroids, particle_coords_cell] = get_particle_centroids(P, particle_list)
% get_particle_centroids computes particle centroids and stores coordinates
% Inputs:
%   P: 3D matrix representing the image with labeled particles
%   particle_list: Nx4 matrix, where the first column contains particle labels,
%               and columns 2-4 contain x, y, z coordinates.
% Outputs:
%   centroids: Matrix with particle ID, count, and centroid coordinates
%   particle_coords_cell: Cell array storing coordinates for each particle

% Get unique particle values and remove zeros
M_unique = unique(P);
M_unique(M_unique == 0) = []; % Remove background

% Preallocate storage
row_c = numel(M_unique);
centroids = zeros(row_c, 5); % Particle ID, count, centroid (y,x,z)
particle_coords_cell = cell(row_c, 1); % Store particle coordinates

% Compute centroids and store coordinates
for c = 1:row_c
    part = M_unique(c);
    idx = particle_list(:, 1) == part;
    part_coords = particle_list(idx, 2:4);
    part_count = sum(idx);

    % Calculate the centroids
    part_y = sum(part_coords(:, 1));
    part_x = sum(part_coords(:, 2));
    part_z = sum(part_coords(:, 3));

    % Store the centroids and count
    centroids(c,1)= part;
    centroids(c,2)= part_count;
    centroids(c,3)= part_y/part_count;
    centroids(c,4)= part_x/part_count;
    centroids(c,5)= part_z/part_count;

    if part_count > 0
        % Store coordinates in the cell array
        particle_coords_cell{c} = part_coords;
    end
end


end
