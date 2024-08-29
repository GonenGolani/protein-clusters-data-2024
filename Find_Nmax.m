function [N_max,index_last] = Find_Nmax(x_max,Nvector,Intensity,cut_off,R_star_index)
%find the last N in the pick at which  Intensity=0
    i=R_star_index;
    x_vec=Nvector./Nvector(R_star_index);
    [~, index_last]=min(abs(x_vec-x_max));
    while i<index_last
        if Intensity(i)<=cut_off
            if Intensity(i)>0
                N_max=Nvector(i);
                index_last=i;
            end
            if Intensity(i)==0
                N_max=Nvector(i-1);
                index_last=i-1;
            end            
            return
        end  
        i=i+1;
    end
    N_max=Nvector(index_last);
end