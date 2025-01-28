function [points, faces] = read_mesh_generic(filename)


    fid = fopen(filename, 'r');
    
    % get the first line
    tline = fgetl(fid);
    tline = strsplit(tline);
        
    
    num_points = str2num(cell2mat(tline(2)));
    num_faces = str2num(cell2mat(tline(4)));
        
   
    points = [];
    for i = 1:num_points
        tline = fgetl(fid);
        tline = strsplit(tline);
        points = [points; str2num(cell2mat(tline(2))) str2num(cell2mat(tline(3))) str2num(cell2mat(tline(4)))];
    end
    
    
    faces = {};
    for i = 1:num_faces
        
        % read this line of face nodes
        tline = fgetl(fid);
        tline = strsplit(tline, {' ','//'}, 'CollapseDelimiters', true);                
        nums = [];
        for j = 1:size(tline,2)            
            this_num = str2num(cell2mat(tline(j)));            
            if min(size(this_num)) > 0                
                nums = [nums this_num];                
            end            
        end
                
        % store this face
        faces{i} = nums;                                         
        
    end
    
    % close the file
    fclose(fid);
    
  
    
    
    
    
    
    

end

