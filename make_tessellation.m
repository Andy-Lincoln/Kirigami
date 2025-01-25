function [pointsD, facesD, points0] = make_tessellation(width, height)
    
    %% we first generate tessellation of initial state
    num_points0 = (1+2*width)*(1+2*height);
    points0 = zeros(num_points0, 2);
    for i = 1:num_points0
        row = ceil(1/(1+2*width));
        col = rem(i,1+2*width);
        if col == 0
            col = 1+2*width;
        end
        points0(i, :) = [col, row];
    end
    % move the pattern such that its center locates at origin
    points0 = points0 - repmat(sum(points0,1)/size(points0,1), size(points0,1),1);
    
    
    num_face0 = 4*width*height;
    face0 = zeros(num_face0, 4);

    for i =1:num_face0
        row = ceil(i/2*width);
        lb_idx = i+row-1;% we can use index of face to determine the index of left bottom vertex
        face0(i, :) = [lb_idx, lb_idx+1, lb_idx+2+2*width, lb_idx+1+2*width];
    end

    %% construct tessellation about deployed state based on the alternate rotation

end