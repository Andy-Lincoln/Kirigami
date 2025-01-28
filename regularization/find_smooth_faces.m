function same_face_adjs = find_smooth_faces(width, height)

    same_face_adjs = zeros(4*(width-1)*(height-1), 2);
    count = 0;
    for i = 1:4*(height-1)*width
        if rem(i+1, 2*width)<=1
            same_face_adjs(count+1, :) = [i, i+4*width];
            count = count+1;
        else
            same_face_adjs(count+1, :) = [i, i+4*width];
            same_face_adjs(count+2, :) = [i, i+2];
            count = count+2;
        end
    end

    for i = 4*(height-1)*width+1: (2*height-1)*2*width-2
        same_face_adjs(count+1, :) = [i, i+2];
        count = count+1;
    end

    for i= (2*height-1)*2*width+1: 4*height*width-2
        same_face_adjs(count+1, :) = [i, i+2];
        count = count+1;
    end

end