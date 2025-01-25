% This is main function

clc, clear

width = 3;
height = 3;
num_points0 = (1+2*width)*(1+2*height);
    points0 = zeros(num_points0, 2);
    for i = 1:num_points0
        row = ceil(i/(1+2*width));
        col = rem(i,1+2*width);
        if col == 0
            col = 1+2*width;
        end
        points0(i, :) = [col, row];
    end
    % move the pattern such that its center locates at origin
    points0 = points0 - repmat(sum(points0,1)/size(points0,1), size(points0,1),1);
    
    
    num_faces = 4*width*height; % the number of faces remains same after deploying
    faces0 = zeros(num_faces, 4);

    for i =1:num_faces
        row = ceil(i/(2*width));
        lb_idx = i+row-1;% we can use index of face to determine the index of left bottom vertex
        faces0(i, :) = [lb_idx, lb_idx+1, lb_idx+2+2*width, lb_idx+1+2*width];
    end
    figure(1)
    clf
    hold on
    axis equal
    axis off
    plot_faces_generic(points0, {faces0}, 1);

    %% construct tessellation about deployed state based on the alternate rotation
    num_pointsD = (8*width+2)*height+2*width;
    pointsD = zeros(num_pointsD, 2);
    Dto0 = zeros(num_pointsD, 1);
    
    facesD = zeros(num_faces, 4);

    pointsD(1:2*width, :) = [1:2:4*width-1; zeros(1, 2*width)]';
    
    % we first construct the faces in the first row
    rotation = repmat([0; 1], width, 1);
    for i=1:2*width
        added_points_ind = 2*width+(2*i+1):-1:2*width+(2*i-1);
        facesD(i,:) = [i, added_points_ind];
        center = [2*i-1, 1];
        pointsD(added_points_ind,:) = [2*i, 1; 2*i-1, 2; 2*i-2, 1];
        Dto0(facesD(i,:)) = circshift(faces0(i, :), rotation(i));
    end
    
    for i = 2:2*height
        


    % figure(2)
    % clf
    % hold on
    % axis equal
    % axis off
    % plot_faces_generic(pointsD(1:19, :), {facesD(1:6, :)}, 2);  
    % 


    % rotation = repmat([1; 0], width, 1); % for constructing map between points0 and pointsD----Dto0

    % facesD(1,:) = [1, 2, 3, 4];
    % pointsD(1:4, :) = [0 1; 1 0; 2 1; 1 2;];
    % Dto0(1:4) = circshift(faces0(1,:), 1);
    % count = 4; % count the number of pointsD constructed
    % 
    % % we first construct the first row of deployed state
    % for i=2:2*width
    %     %center = [2*i-1, 1];
    %     added_vertices = [2*i-1 0; 2*i 1; 2*i-1 2];
    %     pointsD(count+1: count+3, :) = added_vertices;
    %     facesD(i, :) = [3*(i-1) count+1:count+3];
    %     count = count + 3;
    %     Dto0(facesD(i, :)) = circshift(faces0(i, :), rotation(i));
    % end
    % 
    % for i = 2*width+1:num_faces
    %     row = ceil(i/(2*width));
    %     col = rem(i, 2*width);
    %     if col == 0
    %         col = 2*width;
    %     end
    %     center = [2*col-1, 2*row - 1];
    %     if rem(i, 2*width) == 1
    %         added_vertices = [2*col 2*row-1; 2*col-1 2*row; 2*col-2 2*row-1];
    %         pointsD(count+1: count+3, :) = added_vertices;
    %         facesD(i, :) = [lower_face(3) count+1:count+3];
    % 
    % 
    % figure(2)
    % clf
    % hold on
    % axis equal
    % axis off
    % plot_faces_generic(pointsD(1:19, :), {facesD(1:6, :)}, 2);   
