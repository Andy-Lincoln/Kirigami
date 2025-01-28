function [pointsD, face_setsD, edgesD, edge_pairsD, anglesD, int_ringsD, bdy_ringsD, Dto0, points0] = quad_tessellation(width, height, plot_toggle)

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
    
    num_faces = 4*width*height; % the number of faces remains same after deploying
    faces0 = zeros(num_faces, 4);
    
    for i =1:num_faces
        row = ceil(i/(2*width));
        lb_idx = i+row-1;% we can use index of face to determine the index of left bottom vertex
        faces0(i, :) = [lb_idx, lb_idx+1, lb_idx+2+2*width, lb_idx+1+2*width];
    end

    %% construct tessellation about deployed state based on the alternate rotation
    num_pointsD = (8*width+2)*height+2*width;
    pointsD = zeros(num_pointsD, 2);
    Dto0 = zeros(num_pointsD, 1);
    
    facesD = zeros(num_faces, 4);
    
    pointsD(1:2*width, :) = [1:2:4*width-1; zeros(1, 2*width)]';
    
    % we first construct the faces in the first row
    rotation = repmat([-1; 0], width, 1);
    for i=1:2*width
        added_points_ind = 2*width+(2*i+1):-1:2*width+(2*i-1);
        facesD(i,:) = [i, added_points_ind];
        %center = [2*i-1, 1];
        pointsD(added_points_ind,:) = [2*i, 1; 2*i-1, 2; 2*i-2, 1];
        Dto0(facesD(i,:)) = circshift(faces0(i, :), rotation(i));
    end
    count = 2*width; % count the number of facesD constructed
    
    num_M = 2*width + (2*width+1);
    start = 2*width+1;
    for i = 1:2*height-1 % add 2*width+(2*width+1) points for each "M"
        rotation = -1*ones(2*width, 1)-rotation;
        hill_points = facesD(1+(i-1)*2*width:i*2*width, 3);
        added_points_ind = start+i*num_M :1: start+(i+1)*num_M-1;
        pointsD(added_points_ind, :) = pointsD(added_points_ind-num_M, :)+[zeros(num_M, 1) 2*ones(num_M, 1)];
        for j=1:2*width % add 2*width faces to facesD
            count = count+1;
            facesD(count, :) = [hill_points(j), hill_points(j)+num_M+1:-1:hill_points(j)+num_M-1];
            Dto0(facesD(count,:)) = circshift(faces0(count, :), rotation(j));
        end
    end
    
    % move the patterns such that its center locates at origin
    points0 = points0 - repmat(sum(points0,1)/size(points0,1), size(points0,1),1);
    pointsD = pointsD - repmat(sum(pointsD, 1)/size(pointsD, 1), size(pointsD, 1), 1);
    % rescale
    rescale = 1/max(max(abs(pointsD)));
    pointsD = pointsD*rescale;
    points0 = points0*rescale*sqrt(2); % sqrt(2) is the ratio between initial edge (=1) and deployed edge (=sqrt(2))
    face_setsD = {facesD};
    
    %% construct edgesD, edge_pairsD, anglesD, int_ringsD, bdy_ringsD
    edgesD = zeros(4*num_faces, 2);
    anglesD = zeros(4*num_faces, 3); % angle = [A, B, C] representing angle with vertex A and arms/legs AB, AC
    for i=1:num_faces
        face = facesD(i, :);
        edgesD(4*(i-1)+1:4*i, :) = [face' circshift(face, -1)'];
        anglesD(4*(i-1)+1:4*i, :) = [face' circshift(face, -1)' circshift(face, 1)'];
    end
    
    % determine edge pair
    edge_pairsD = make_edge_pairs(edgesD, Dto0);
    
    % determine interior_rings: four angles in anglesD sum up to 2*pi
    % and boundary_rings: one angle = pi/2 or two angles sum up to pi
    
    % we first rearrange anglesD according to vertex of its corresponding initial
    % angle, which facilates finding int_rings and bdy_rings greatly thanks to
    % the order of points0 we constructed!
    [~, ind] = sortrows(Dto0(anglesD));
    anglesD = anglesD(ind,:);
    
    int_ringsD = zeros((2*width-1)*(2*height-1), 2);
    bdy_ringsD = zeros(4*(width+height), 3);
    bdy_ringsD(1, :) = [1 1 pi/2];
    count = 1;
    for i=2:2*width
        bdy_ringsD(i, :) = [count+1 count+2 pi];
        count = count+2;
    end
    bdy_ringsD(2*width+1, :) = [count+1 count+1 pi/2];
    count = count+1;
    tmp1 = 2*width+1;
    tmp2 = 0;
    for i=2*width+2: 2*height*(2*width+1)
        if rem(i, 2*width+1) <= 1
            tmp1 = tmp1+1;
            bdy_ringsD(tmp1,:) = [count+1 count+2 pi];
            count = count+2;
        else
            tmp2 = tmp2+1;
            int_ringsD(tmp2, :) = [count+1, count+4];
            count = count + 4;
        end
    end
    bdy_ringsD(tmp1+1, :) = [count+1 count+1 pi/2];
    count = count+1;
    for i=2:2*width
        bdy_ringsD(tmp1+i, :) = [count+1 count+2 pi];
        count = count+2;
    end
    bdy_ringsD(end, :) = [4*num_faces, 4*num_faces, pi/2];
   
    %% draw edge_pairsD, ringsD for illustration
    figure(1)
    clf
    hold on
    axis equal
    axis off
    plot_faces_generic(points0, {faces0}, 1);
    
    figure(2)
    clf
    hold on
    axis equal
    axis off
    plot_faces_generic(pointsD, face_setsD, 2);
    
    if plot_toggle
        figure(2)
        hold on
        for i = 1:size(edge_pairsD,1)
            ptA = (pointsD(edgesD(edge_pairsD(i,1),1),:) + pointsD(edgesD(edge_pairsD(i,1),2),:))/2;
            ptB = (pointsD(edgesD(edge_pairsD(i,2),1),:) + pointsD(edgesD(edge_pairsD(i,2),2),:))/2;
            plot([ptA(1) ptB(1)], [ptA(2) ptB(2)], '--r', 'LineWidth',2)
        end
    end

end