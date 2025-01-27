function pairs = make_edge_pairs(edges, Dto0)

    num_edges = size(edges,1);
    pairs = [];            
    paired = zeros(num_edges,1);
    for i = 1:num_edges-1       
        if paired(i) == 1
            continue
        end
        this_edge = edges(i,:);                 
        for j = i+1:num_edges
            if paired(j) == 1
                continue
            end
            other_edge = edges(j,:);     
            % one cannot use Dto0(this_edge)==Dto0(other_edge) to judge 
            % for example [2 7] [7 2] are indeed same edge.
            if sum(ismember(Dto0(this_edge), Dto0(other_edge))) == 2                                             
                pairs = [pairs; i j];
                paired(i) = 1;
                paired(j) = 1;
                break                
            end            
        end
    end

end

