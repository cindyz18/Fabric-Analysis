function B = map_porosity_to_matrix(porosity, a1, a2, a3, l)
% map_porosity_to_matrix maps porosity values to a 4D matrix efficiently
%
% Inputs:
%   porosity: Column vector of porosity values
%   a1, a2, a3: Dimensions of the target 4D matrix
%   l: Adjustment factor for matrix dimensions
%
% Output:
%   B: 4D matrix with mapped porosity values
%
% Description:
%   This function normalizes porosity values into a 0-255 range, then efficiently maps them 
%   to a 4D matrix without using explicit loops.

% Normalize porosity values into a 0-255 range
min_val = min(porosity);
max_val = max(porosity);
bin_range = max_val - min_val;
A = [porosity, round(((porosity - min_val) * (255-0) / (bin_range))+0, 0)];

disp("Completed mapping, reshaping matrix");

% % Reshape directly into a 4D matrix
B = reshape(A(:,2), [(a3-l), (a2-l), (a1-l), 1]); % Flipped dimensions as it will be rearranged
B = permute(B, [3, 2, 1]);

% Convert to uint16
B = uint16(B);

end
