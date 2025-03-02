function porosity = calculate_porosity_values(M1, L)
% calculate_stored_values computes porosity values using the movsum method
% Inputs:
%   M1: 3D matrix from which porosity is calculated
%   L: Size of the region of interest (ROI) for movsum operation
% Outputs:
%   stored: Column vector of porosity values

[a1, a2, a3] = size(M1);
moving_sum_i = movsum(M1, L, 1, 'Endpoints','discard');
moving_sum_j = movsum(moving_sum_i, L, 2, 'Endpoints','discard');
moving_sum_k = movsum(moving_sum_j, L, 3, 'Endpoints','discard');
porosity = permute(moving_sum_k(1:a1-L,1:a2-L,1:a3-L),[3 2 1]);
porosity = porosity(:);

end