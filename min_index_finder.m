function [index_out] = min_index_finder(x_vector,index_at_R_star,x_min)
    [~,index_x_min]=min(abs(x_vector-x_min));
    if index_x_min>index_at_R_star
        index_out=index_x_min;
    else
        index_out=index_at_R_star;
    end

end