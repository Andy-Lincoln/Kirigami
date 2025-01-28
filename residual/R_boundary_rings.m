function R = R_boundary_rings(points, angles, bdy_rings)

    % bdy_rings: n x 3, [a,b,c] anglesD from a to b sum up to c (pi or pi/2)
    % dimensions    
    num_rings = size(bdy_rings, 1);
    
    % initialize R
    R = zeros(num_rings,1);
    
    % compute and store each set of boundary ring residuals 
    for i = 1:num_rings
        
        % loop through angles in this ring
        for j = bdy_rings(i,1):bdy_rings(i,2)
            
            % compute this angle            
            this_angle = atan2_angle(points(angles(j,1),:), ...
                                     points(angles(j,2),:), ...
                                     points(angles(j,3),:));          
            
            % store this angle
            R(i) = R(i) + this_angle;
            
        end
        
        % store this angle
        R(i) = R(i) - bdy_rings(i, 3);
        
    end    
        
end

