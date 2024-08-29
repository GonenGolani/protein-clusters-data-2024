function [max_x,Max_Value] = Find_R_Star(X_value,y_value,max_index)
xq=linspace(X_value(max_index-4),X_value(max_index+4),100);
s = spline(X_value,y_value,xq);

[Max_Value,index]=max(s);
max_x=xq(index);
end