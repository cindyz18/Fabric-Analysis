function [particle_list, contact_list, M] = get_contact_list(P, C)
    % calculate_stored_values computes porosity values using the movsum method
    % Inputs:
    %   L: 3D matrix that defines the labelled particles
    %   C: 3D matrix that defines the contacts
    % Outputs:
    %   particle_list: Matrix of the particles
    %   contact_list: Matrix of the contacts
    %   M: manipulated matrix

    InitalConVoxVal = max(C ,[], "all" ); % Gets the maximum voxel value from the contacts
    ConVoxVal = max(P ,[], "all" ) + 1; % Gets the maximum voxel value from the labelled particles
    disp(['InitalConVoxVal: ', num2str(InitalConVoxVal), ', ConVoxVal: ', num2str(ConVoxVal)]);   
    
    % Changing all values of contact voxels to a constant value of ConVoxVal -
    % This is done as the contacts in Avizo cannot be assigned value above 255
    M = P;
    M(C==InitalConVoxVal) = ConVoxVal;
    
    [rows, columns, slices] = size(M);
    fprintf('rows: %d, columns: %d, slices: %d\n',rows,columns,slices)
    
    %% Optimised code
    tic;
    disp('=========================================================')
    % Flattened matrix
    M_list= zeros(numel(M), 4);                            %numel is the no.of elements in the matrix 
    
    % Use vectorization to avoid nested loops
    [X, Y, Z] = ndgrid(1:size(M,1), 1:size(M,2), 1:size(M,3));
    M_flat = M(:);
    
    % Fill M_list using the precomputed indices and M values
    M_list(:, 1) = M_flat;
    M_list(:, 2) = X(:);
    M_list(:, 3) = Y(:);
    M_list(:, 4) = Z(:); 
    
    disp("value,row,column,slice")
    
    first_loop_time = toc;
    disp(['Time taken by first_loop: ', num2str(first_loop_time), ' seconds']);
  
    tic;
    % Logical indexing to find indices where the condition is met
    idx = M_list(:, 1) > ConVoxVal;
    % Update the values in M_list that meet the condition
    M_list(idx, 1) = ConVoxVal;
    
    second_loop_time = toc;
    disp(['Time taken by second loop: ', num2str(second_loop_time), ' seconds']);
    
    
    % Removing all the coordinates of '0' and 'ConVoxVal' so that only particles are visible
    % Define the conditions for just labelled particles, not including the contacts
    condition = (M_list(:, 1) > 0) & (M_list(:, 1) < ConVoxVal);
    % Use logical indexing to select rows that satisfy the condition
    particle_list = M_list(condition, :);
    
    % Define the condition for the contacts
    contact_condition = M_list(:, 1) == ConVoxVal;
    
    % Use logical indexing to select rows that satisfy the condition
    contact_list = M_list(contact_condition, :);

end