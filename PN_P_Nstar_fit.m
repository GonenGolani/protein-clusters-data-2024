function [y] = PN_P_Nstar_fit(x,C3,C2,C1,data_type)
    if strcmp(data_type,'Number')
        y=exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
    end
    if strcmp(data_type,'Volume')
        y=x.^3.*exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
    end
    if strcmp(data_type,'Intensity')
        y=x.^6.*exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
    end 
end