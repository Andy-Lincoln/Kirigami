function [rings, anglesD] = make_dev_rings_generic(anglesD, Dto0)
% returns a sorted anglesD and intervals constituting rings on anglesD
    
    % map the deployed angles to initial indices
    angles0 = Dto0(anglesD);
    
    % sort both by central node in angles
    [angles0, ind] = sortrows(angles0);    
    anglesD = anglesD(ind,:);      
    
    % find angles belonging to the same central node
    ring_breaks = find([true; logical(diff(angles0(:,1)))]);
    ring_breaks = [ring_breaks; size(angles0,1)+1];
    
    % determine if outer nodes in angles at each central node form a ring
    rings = [];
    for i = 1:length(ring_breaks)-1
       
        % get the outer nodes
        this_ring = angles0(ring_breaks(i):ring_breaks(i+1)-1, [2 3]);
        outer_nodes = sort(reshape(this_ring, 2*size(this_ring, 1),1));
        
        % check whether each outer node appears exactly twice
        if checkDuplicate(outer_nodes)            
            rings = [rings; ring_breaks(i) ring_breaks(i+1)-1];
        end
        
    end
        
end


function dup = checkDuplicate(a)
    dup = 1;
    for i=1:2:length(a)-1
        if a(i) ~= a(i+1) % we don't check a(i)==a(i+2) because it's impossible in our case, which means this function is not versatile!
            dup = 0;
            break
        end
    end
end

