clear
clc
tic;
% file_name = 'UG_6pc_sand_binary_v1';
% file_path = 'C:\Users\zcin1\OneDrive\Documents\Uni\2024\FYP\final samples\MATLAB\UG\B\';
file_name_matrix = {'UG_0pc_sand_binary_v1' 1};

for i = 1:size(file_name_matrix,1)
    tic;
    file_name = file_name_matrix{i,1};
    display(['file name ', file_name]);
    load(strcat(file_name, '.mat'));   %Evaluate the first 3 lines - then copy the generated matrix name into line 4.
    eval(strcat('M0=Avizo_',file_name,'_mat;')) %Do not remove or change the first 6 and last 4 characters - "Avizo_" and "_mat"
    M0=squeeze(M0);
    % M0(:,:,1)=[1 2 3 123 321; 4 5 6 456 654; 7 8 9 789 987; 1 2 3 4 5; 6 7 8 9 10];
    % M0(:,:,2)=[10 11 12 111 112; 13 14 15 345 543; 16 17 18 678 876; 5 6 7 8 9; 6 7 8 9 10];
    % M0(:,:,3)=[19 20 21 901 109; 22 23 24 234 432; 25 26 27 567 765; 9 10 11 12 13; 6 7 8 9 10];
    % M0(:,:,4)=[28 29 30 890 908; 4 14 6 446 644; 7 8 9 789 987; 1 2 3 4 5; 6 7 8 9 10];
    % M0(:,:,5)=[10 11 12 112 211; 13 14 15 345 543; 16 17 18 678 876; 5 6 7 8 9; 6 7 8 9 10];
    % M0(:,:,6)=[19 20 21 231 132; 22 23 24 234 432; 25 26 27 567 765; 9 10 11 12 13; 6 7 8 9 10];
    
    %%My inputs
    REV = 3:10:250; %Increase the matrix to suit
    pos = [350 350 505; 350 650 505; 650 350 505; 650 650 505];
    [y1,y2] = size(REV);
    [x1,x2] = size(pos);
    M2 = zeros(x1*y2,5);        %matrix to store the porosity values for different positions but for one REV
    counter=1;
    for jj=1:y2
        len = REV(1,jj)-1    
        for ii=1:x1                                                          
            M1 = M0(pos(ii,1)-(len/2):pos(ii,1)+(len/2),pos(ii,2)-(len/2):pos(ii,2)+(len/2),pos(ii,3)-(len/2):pos(ii,3)+(len/2));        %Slicing the matrix
            [a1, a2, a3] = size(M1);
            L=1;                                                                                         %por_roi size 
            stored=zeros(((a1-L)*(a2-L)*(a3-L)),1);
            c=1;
            inc=1;
            for i=1:inc:(a1-L);
                for j=1:inc:(a2-L);
                    for k=1:inc:(a3-L);
                        por_roi=M1(i:i+(L-1),j:j+(L-1),k:k+(L-1));
                        a=sum(por_roi,1);
                        f=sum(a);
                        n=sum(f);
                        ratio=n/(L*L*L);
                        stored(c,1)=ratio;
                        c=c+1;                  
                    end
                end
            end
            stored;
            porosity_avg=sum(stored)/((a1-L)*(a2-L)*(a3-L));
            M2(counter,:)= [REV(jj),pos(ii,1),pos(ii,2),pos(ii,3),porosity_avg];             %REV,Position,Porosity
            counter=counter+1;
            
        end 
    end
    col_header={'REV size','Loc_x','Loc_y','Loc_z','Porosity'}; 
    xlswrite('REVfinder.xlsx',M2,file_name,'B2');     %Write data
    xlswrite('REVfinder.xlsx',col_header,file_name,'B1');     %Write column header
    time_taken = toc;
    disp(['Time taken: ', num2str(time_taken), ' seconds']);
end