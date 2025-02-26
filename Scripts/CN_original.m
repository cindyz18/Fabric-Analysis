clear;clc;clear all;close all;
csv_file_path = '/fs03/uc40/Cindy/final_samples/CN/CN_csv/';
mat_file_path = '/fs03/uc40/Cindy/final_samples/CN/mat_files/';

addpath('./../Functions','./../Data')
file_names = {'UG_0pc_sand_contact_binary_v1' 'UG_0pc_sand_separated_labelled_v1'};
%     ;...
% 'UG_3pc_sand_contact_binary_v1' 'UG_3pc_sand_separated_labelled_v1';...
% 'UG_6pc_sand_contact_binary_v1' 'UG_6pc_sand_separated_labelled_v2';...
% 'UG_9pc_sand_contact_binary_v1' 'UG_9pc_sand_separated_labelled_v1';...
% 'UG_15pc_sand_contact_binary_v1' 'UG_15pc_sand_separated_labelled_v1';...
% 'GP_15pc_sand_contact_binary_v1' 'GP_15pc_sand_binary_separated_labelled_v1';...
% 'GP_9pc_sand_contact_binary_v1' 'GP_9pc_sand_binary_separated_labelled_v1';...
% 'GP_6pc_sand_contact_binary_manual_v1' 'GP_6pc_sand_binary_manual_separated_labelled_v1';...
% 'GP_3pc_sand_contact_binary_v1' 'GP_3pc_sand_binary_separated_labelled_v1';...
% 'GP_0p_sand_contact_binary_v1' 'GP_0p_sand_binary_separated_labelled_v1'};

for file_no = 1:size(file_names,1)
    contact_file_name = file_names{file_no,1};
    sand_file_name = file_names{file_no,2};
    
    load(strcat(contact_file_name, '.mat')); %Evaluate the first 3 lines - then copy the generated matrix name into line 4.
    eval(strcat('contacts=Avizo_',contact_file_name,'_mat;'))  %Do not remove or change the first 6 and last 4 characters - "Avizo_" and "_mat"
    contacts=squeeze(contacts);
    
    load(strcat(sand_file_name, '.mat'));
    eval(strcat('labelled=Avizo_',sand_file_name,'_mat;'))
    labelled=squeeze(labelled);

    contacts = contacts(201:400,201:400,201:400);
    labelled = labelled(201:400,201:400,201:400);
    
    InitalConVoxVal = max(contacts ,[], "all" );
    ConVoxVal = max(labelled ,[], "all" )+1;
    disp(['InitalConVoxVal: ', num2str(InitalConVoxVal), ', ConVoxVal: ', num2str(ConVoxVal)]);
    
    % contacts = contacts(201:300,201:300,201:300);
    % labelled = labelled(201:300,201:300,201:300);  

     
    
    M0 = labelled;
    M0(contacts==InitalConVoxVal) = ConVoxVal;
    
    [rows, columns, slices] = size(M0);
    fprintf('rows: %d, columns: %d, slices: %d\n',rows,columns,slices)
    
    %% Optimised code
    tic;
    disp('=========================================================')
    M_list_opt= zeros(numel(M0), 4);                            %numel is the no.of elements in the matrix 
    counter = 1;
    
    % Use vectorization to avoid nested loops
    [X, Y, Z] = ndgrid(1:size(M0,1), 1:size(M0,2), 1:size(M0,3));
    M0_flat = M0(:);
    
    % Fill M_list_opt using the precomputed indices and M0 values
    M_list_opt(:, 1) = M0_flat;
    M_list_opt(:, 2) = X(:);
    M_list_opt(:, 3) = Y(:);
    M_list_opt(:, 4) = Z(:); 
    
    [m1,m2]=size(M_list_opt);
    disp("value,row,column,slice")
    
    first_loop_time = toc;
    disp(['Time taken by first_loop: ', num2str(first_loop_time), ' seconds']);
    
    %%%Changing all values of contact voxels to a constant value of ConVoxVal -
    %%%This is done as the contacts in Avizo cannot be assigned value above 255
    tic;
    
    % Logical indexing to find indices where the condition is met
    idx = M_list_opt(:, 1) > ConVoxVal;
    % Update the values in M_list_opt that meet the condition
    M_list_opt(idx, 1) = ConVoxVal;
    
    second_loop_time = toc;
    disp(['Time taken by second loop: ', num2str(second_loop_time), ' seconds']);
    
    %%%%%%%%%%%%%%%%%%%
    %%% VISUAL PART %%%
    %%%%%%%%%%%%%%%%%%%
    
    %%%Removing all the coordinates of '0' and 'ConVoxVal' so that only particles are visible
    
    tic;
    
    % Define the conditions
    condition = (M_list_opt(:, 1) > 0) & (M_list_opt(:, 1) < ConVoxVal);
    
    % Use logical indexing to select rows that satisfy the condition
    M_refinedlist_opt = M_list_opt(condition, :);
    
    third_loop_time = toc;
    disp(['Time taken by third loop: ', num2str(third_loop_time), ' seconds']);
    % fprintf('Is equal\nM_refinedlist,M_refinedlist_opt: %d \n', isequal(M_refinedlist,M_refinedlist_opt))
    
    
    tic;
    
    % Define the condition
    contact_condition = M_list_opt(:, 1) == ConVoxVal;
    
    % Use logical indexing to select rows that satisfy the condition
    M_contactlist_opt = M_list_opt(contact_condition, :);
    
    
    % figure(file_no);
    % plot3(M_contactlist_opt(:,2),M_contactlist_opt(:,3),M_contactlist_opt(:,4),'.','MarkerSize',20,'Color','yellow');       
    % hold on
    
    fourth_loop = toc;
    disp(['Time taken by fourth loop: ', num2str(fourth_loop), ' seconds']);
    % fprintf('Is equal\nM_contactlist,M_contactlist_opt: %d \n', isequal(M_contactlist,M_contactlist_opt))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Getting the 26 neighbourhood values of a contact voxel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    
    [a1, a2, a3] = size(M0);
    count_M0 = 1;
    % Preallocate for maximum possible size
    M_C_location_opt = zeros(size(M_contactlist_opt,1), 3);
    M_C_neighbours_opt = zeros(size(M_contactlist_opt,1), 26);
    
    for idx = 1:size(M_contactlist_opt,1)
        i = M_contactlist_opt(idx, 2);
        j = M_contactlist_opt(idx, 3);
        k = M_contactlist_opt(idx, 4);
        % Ensure we are not at the boundary to avoid out-of-bounds error
        if i > 1 && i < a1 && j > 1 && j < a2 && k > 1 && k < a3
            % Extract the 3x3x3 neighborhood - getting the 27 values - 26 neighbours and R_voxel value in a matrix
            neighbors = M0(i-1:i+1, j-1:j+1, k-1:k+1);
            
            % Reshape and remove the central voxel value
            % Reshaping the matrix into a single row array
            Contact_neighbors = reshape(neighbors, 1, []);
            Contact_neighbors(14) = [];
            
            % Store the location and neighbors
            M_C_location_opt(count_M0, :) = [i, j, k];
            M_C_neighbours_opt(count_M0, :) = Contact_neighbors;
            
            count_M0 = count_M0 + 1;
        end
    
    end
    
    % Remove unused preallocated rows
    M_C_location_opt = M_C_location_opt(1:count_M0-1, :);
    M_C_neighbours_opt = M_C_neighbours_opt(1:count_M0-1, :);
    
    M_C_location_opt;
    M_C_neighbours_opt;
    
    fifth_loop = toc;
    disp(['Time taken by fifth loop: ', num2str(fifth_loop), ' seconds']);
    % fprintf('Is equal\nM_C_location,M_C_location_opt: %d \n', isequal(M_C_location,M_C_location_opt))
    % fprintf('Is equal\nM_C_neighbours,M_C_neighbours_opt: %d \n', isequal(M_C_neighbours,M_C_neighbours_opt))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% COORDINATION NUMBER %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%Getting the unique values of the rows and storing them
    tic;
    
    % Preallocate the output matrix M3 with zeros
    [a1, a2] = size(M_C_neighbours_opt);
    M3_opt = zeros(a1, 2); 
    
    % Replace zeros with ConVoxVal and get unique values row-wise
    M2_modified_opt = M_C_neighbours_opt;
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
    M3_uniq_opt;
    
    a = M3_uniq_opt;                          %Counting algorithm
    [voxel_value,ia,ic] = unique(a);
    a_counts_opt = accumarray(ic,1);
    value_counts_opt = [voxel_value, a_counts_opt]; 
    
    disp("Coordination number : [Voxel value, Repetitions]")
    % value_counts_opt
    save(strcat(mat_file_path, sand_file_name,'_value_counts_opt_v3.mat'),'value_counts_opt','-v7.3')
    col_header={'Voxel value', 'Repetitions'}; 
    xlswrite(strcat(csv_file_path,sand_file_name,'_CN.xlsx'),value_counts_opt,'CN','B2');     %Write data
    % xlswrite(strcat(sand_file_name,'_CN.xlsx'),col_header,'CN','B1');     %Write column header
    
    sixth_loop = toc;
    disp(['Time taken by sixth loop: ', num2str(sixth_loop), ' seconds']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Getting the Centroids of all the particles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tic;
    % 
    % M_list_opt = M_refinedlist_opt;
    % [m1, ~] = size(M_list_opt);
    % 
    % % Get unique values in the image and remove zeros
    % M_uniQ_opt = unique(M0);
    % M_uniQ_opt(M_uniQ_opt == 0) = [];
    % 
    % 
    % [row_c, ~] = size(M_uniQ_opt);
    % 
    % % Preallocate for centroids
    % M_centroids_opt = zeros(row_c, 5);
    % 
    % for c = 1:row_c
    %     part = M_uniQ_opt(c);
    % 
    %     % Find indices where the particle matches
    %     idx = M_list_opt(:, 1) == part;
    % 
    %     % Get the coordinates and count
    %     part_coords = M_list_opt(idx, 2:4);
    %     part_count = sum(idx);
    % 
    %     % Calculate the centroids
    %     part_y = sum(part_coords(:, 1));
    %     part_x = sum(part_coords(:, 2));
    %     part_z = sum(part_coords(:, 3));
    % 
    %     % Store the centroids and count
    %     M_centroids_opt(c,1)= part;
    %     M_centroids_opt(c,2)= part_count;
    %     M_centroids_opt(c,3)= part_y/part_count;
    %     M_centroids_opt(c,4)= part_x/part_count;
    %     M_centroids_opt(c,5)= part_z/part_count;
    % 
    %     % % Plot individual particles
    %     % if part_count > 0
    %     %     k = boundary(part_coords(:, 1), part_coords(:, 2), part_coords(:, 3));
    %     %     colorbar
    %     %     hold on;
    %     %     trisurf(k, part_coords(:, 1), part_coords(:, 2), part_coords(:, 3), 'Facecolor', 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none');
    %     % end
    % end
    % 
    % M_centroids_opt; % Particle, Voxelcount, cg_y, cg_x, cg_z
    % save(strcat(sand_file_name,'_M_centroids_opt_v3.mat'),'M_centroids_opt','-v7.3')
    % col_header={'Particle', 'Voxelcount', 'cg_y', 'cg_x', 'cg_z'}; 
    % xlswrite(strcat(sand_file_name,'_CN.xlsx'),M_centroids_opt,'centroid','B2');     %Write data
    % xlswrite(strcat(sand_file_name,'_CN.xlsx'),col_header,'centroid','B1');     %Write column header
    % 
    % % % Plot the centroids of particles as points
    % % plot3(M_centroids_opt(:, 3), M_centroids_opt(:, 4), M_centroids_opt(:, 5), '.', 'MarkerSize', 20, 'Color', 'green');
    % 
    % seventh_loop = toc;
    % 
    % disp(['Time taken by seventh_loop: ', num2str(seventh_loop), ' seconds']);

end
