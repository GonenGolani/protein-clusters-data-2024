function [C3_out,C2_out,C1_out,Cm1_out , chi2_min,Error_RMS] = LSM_C1C2C3_no_constreint(x_all,y_all,Index_first,index_last,C_limit,dC,fit_type,data_type,with_Neg_C1)
x=x_all(Index_first:index_last);
y=y_all(Index_first:index_last);
N_point=length(x);

i=1;
C2=0;
while C2<=C_limit
if with_Neg_C1
    C1=-C_limit;
else
    C1=0;
end
    while C1<C_limit
        C3=0;
        while C3<C_limit
            C1_vec(i)=C1;
            C2_vec(i)=C2;
            C3_vec(i)=C3;
            if strcmp(data_type,'Number')
                y_model=exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
            end
            if strcmp(data_type,'Volume')
                y_model=x.^3.*exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
            end
            if strcmp(data_type,'Intensity')
                y_model=x.^6.*exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
            end 
            if strcmp(fit_type,'exp')
                chi2(i)=sum((y-y_model).^2);
                error(i)=(sum(abs(y-y_model)));
            end
            if strcmp(fit_type,'log')

                chi2(i)=sum((log(y)-log(y_model)).^2);
                error(i)=(sum(abs(log(y)-log(y_model))));

            end
        
            i=i+1;
            C3=C3+dC*10;
        end
        C1=C1+dC*10;
    end
   C2=C2+dC*10;
end
[chi2_min,min_index]=min(chi2);
C3_out=C3_vec(min_index);
C2_out=C2_vec(min_index);
C1_out=C1_vec(min_index);

%% round 2

i=1;
C2=max(C2_out-dC*10,0);
while C2<=C2_out+dC*1
if with_Neg_C1
    C1=C1_out-dC*10;
else
    C1=max(C1_out-dC*10,0);
end
    while C1<C1_out+dC*10
        C3=max(C3_out-dC*10,0);
        while C3<C3_out+dC*10
            C1_vec(i)=C1;
            C2_vec(i)=C2;
            C3_vec(i)=C3;
            if strcmp(data_type,'Number')
                y_model=exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
            end
            if strcmp(data_type,'Volume')
                y_model=x.^3.*exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
            end
            if strcmp(data_type,'Intensity')
                y_model=x.^6.*exp(-C1*(x.^1-1)-C2*(x.^2-1)-C3*(x.^3-1));
            end 
            if strcmp(fit_type,'exp')
                chi2(i)=sum((y-y_model).^2);
                error(i)=(sum(abs(y-y_model)));
            end
            if strcmp(fit_type,'log')

                chi2(i)=sum((log(y)-log(y_model)).^2);
                error(i)=(sum(abs(log(y)-log(y_model))));

            end
        
            i=i+1;
            C3=C3+dC;
        end
        C1=C1+dC;
    end
   C2=C2+dC;
end
[chi2_min,min_index]=min(chi2);
C3_out=C3_vec(min_index);
C2_out=C2_vec(min_index);
C1_out=C1_vec(min_index);

Cm1_out=0;
Error_RMS=error(min_index)/(N_point-3);


end





