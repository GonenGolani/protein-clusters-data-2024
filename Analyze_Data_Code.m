%% plot all proteins normelized drop size
clc
clear all
close all
restoredefaultpath;

%% what to plot
FUS_NTA=1;
FUS_DLS=1;
CPEB4=1;
hnRNPA3=1;
EWSR1=0;
TAF15=0;

dC=0.01;

y_min=0.0001;
x_min=1.0;
x_max=5;


KB=1.38*10^-23;
x_vec=linspace(x_min,x_max,1000);

%% colors

Col_025=[228 26 28]/256;	
Col_050=[55 126 184]/256;	
Col_100=[77 175 74]/256;	
Col_200=[152 78 163]/256;
Col_300='k';	
Col_0125=[255 127 0]/256;
Col_600='k';
Col_070=[166 86 40]/256;

%% what to fit

if FUS_NTA
Data_type_FUS_NTA='Number';
data_NTA_FUS=importdata('FUS NTA.txt');

Radius_FUS_NTA=data_NTA_FUS(:,1)/2*10^-9;
FUS_NTA_0125=data_NTA_FUS(:,2)./sum(data_NTA_FUS(:,2)); % 0.25 micro molar
FUS_NTA_025=data_NTA_FUS(:,3)./sum(data_NTA_FUS(:,3)); % 0.50 micro molar
FUS_NTA_050=data_NTA_FUS(:,4)./sum(data_NTA_FUS(:,4)); % 0.70 micro molar
FUS_NTA_100=data_NTA_FUS(:,5)./sum(data_NTA_FUS(:,5)); % 1.00 micro molar
FUS_NTA_200=data_NTA_FUS(:,6)./sum(data_NTA_FUS(:,6)); % 2.00 micro molar


% find radius and distribution width
[~,index_FUS_NTA(1)]=max(FUS_NTA_0125); 
[~,index_FUS_NTA(2)]=max(FUS_NTA_025); 
[~,index_FUS_NTA(3)]=max(FUS_NTA_050);
[~,index_FUS_NTA(4)]=max(FUS_NTA_100); 
[~,index_FUS_NTA(5)]=max(FUS_NTA_200);

[R_star_FUS_NTA(1),FUS_NTA_max_vec(1)] =Find_R_Star(Radius_FUS_NTA,FUS_NTA_0125,index_FUS_NTA(1));
[R_star_FUS_NTA(2),FUS_NTA_max_vec(2)] =Find_R_Star(Radius_FUS_NTA,FUS_NTA_025,index_FUS_NTA(2));
[R_star_FUS_NTA(3),FUS_NTA_max_vec(3)] =Find_R_Star(Radius_FUS_NTA,FUS_NTA_050,index_FUS_NTA(3));
[R_star_FUS_NTA(4),FUS_NTA_max_vec(4)] =Find_R_Star(Radius_FUS_NTA,FUS_NTA_100,index_FUS_NTA(4));
[R_star_FUS_NTA(5),FUS_NTA_max_vec(5)] =Find_R_Star(Radius_FUS_NTA,FUS_NTA_200,index_FUS_NTA(5));


x_FUS_NTA_0125=Radius_FUS_NTA/R_star_FUS_NTA(1);
x_FUS_NTA_025=Radius_FUS_NTA/R_star_FUS_NTA(2);
x_FUS_NTA_050=Radius_FUS_NTA/R_star_FUS_NTA(3);
x_FUS_NTA_100=Radius_FUS_NTA/R_star_FUS_NTA(4);
x_FUS_NTA_200=Radius_FUS_NTA/R_star_FUS_NTA(5);

min_index_FUS_NTA(1)=min_index_finder(x_FUS_NTA_0125,index_FUS_NTA(1),x_min);
min_index_FUS_NTA(2)=min_index_finder(x_FUS_NTA_025,index_FUS_NTA(2),x_min);
min_index_FUS_NTA(3)=min_index_finder(x_FUS_NTA_050,index_FUS_NTA(3),x_min);
min_index_FUS_NTA(4)=min_index_finder(x_FUS_NTA_100,index_FUS_NTA(4),x_min);
min_index_FUS_NTA(5)=min_index_finder(x_FUS_NTA_200,index_FUS_NTA(5),x_min);

y_FUS_NTA_0125=FUS_NTA_0125./FUS_NTA_max_vec(1);
y_FUS_NTA_025=FUS_NTA_025./FUS_NTA_max_vec(2);
y_FUS_NTA_050=FUS_NTA_050./FUS_NTA_max_vec(3);
y_FUS_NTA_100=FUS_NTA_100./FUS_NTA_max_vec(4);
y_FUS_NTA_200=FUS_NTA_200./FUS_NTA_max_vec(5);

[~,FUS_NTA_index_last(1)] = Find_Nmax(x_max,Radius_FUS_NTA,y_FUS_NTA_0125,y_min,index_FUS_NTA(1));
[~,FUS_NTA_index_last(2)] = Find_Nmax(x_max,Radius_FUS_NTA,y_FUS_NTA_025,y_min,index_FUS_NTA(2));
[~,FUS_NTA_index_last(3)] = Find_Nmax(x_max,Radius_FUS_NTA,y_FUS_NTA_050,y_min,index_FUS_NTA(3));
[~,FUS_NTA_index_last(4)] = Find_Nmax(x_max,Radius_FUS_NTA,y_FUS_NTA_100,y_min,index_FUS_NTA(4));
[~,FUS_NTA_index_last(5)] = Find_Nmax(x_max,Radius_FUS_NTA,y_FUS_NTA_200,y_min,index_FUS_NTA(5));

[FUS_NTA_C3(1),FUS_NTA_C2(1),FUS_NTA_C1(1),FUS_NTA_Cm1(1),FUS_NTA_chi2(1),FUS_NTA_Error_RMS(1)] = LSM_C1C2C3_no_constreint(x_FUS_NTA_0125,y_FUS_NTA_0125,min_index_FUS_NTA(1),FUS_NTA_index_last(1),8,dC,'log',Data_type_FUS_NTA,0);
[FUS_NTA_C3(2),FUS_NTA_C2(2),FUS_NTA_C1(2),FUS_NTA_Cm1(2),FUS_NTA_chi2(2),FUS_NTA_Error_RMS(2)] = LSM_C1C2C3_no_constreint(x_FUS_NTA_025,y_FUS_NTA_025,min_index_FUS_NTA(2),FUS_NTA_index_last(2),8,dC,'log',Data_type_FUS_NTA,0);
[FUS_NTA_C3(3),FUS_NTA_C2(3),FUS_NTA_C1(3),FUS_NTA_Cm1(3),FUS_NTA_chi2(3),FUS_NTA_Error_RMS(3)] = LSM_C1C2C3_no_constreint(x_FUS_NTA_050,y_FUS_NTA_050,min_index_FUS_NTA(3),FUS_NTA_index_last(3),8,dC,'log',Data_type_FUS_NTA,0);
[FUS_NTA_C3(4),FUS_NTA_C2(4),FUS_NTA_C1(4),FUS_NTA_Cm1(4),FUS_NTA_chi2(4),FUS_NTA_Error_RMS(4)] = LSM_C1C2C3_no_constreint(x_FUS_NTA_100,y_FUS_NTA_100,min_index_FUS_NTA(4),FUS_NTA_index_last(4),8,dC,'log',Data_type_FUS_NTA,0);
[FUS_NTA_C3(5),FUS_NTA_C2(5),FUS_NTA_C1(5),FUS_NTA_Cm1(5),FUS_NTA_chi2(5),FUS_NTA_Error_RMS(5)] = LSM_C1C2C3_no_constreint(x_FUS_NTA_200,y_FUS_NTA_200,min_index_FUS_NTA(5),FUS_NTA_index_last(5),8,dC,'log',Data_type_FUS_NTA,0);


y_FUS_NTA_0125_fit=PN_P_Nstar_fit(x_vec,FUS_NTA_C3(1),FUS_NTA_C2(1),FUS_NTA_C1(1),Data_type_FUS_NTA);
y_FUS_NTA_025_fit=PN_P_Nstar_fit(x_vec,FUS_NTA_C3(2),FUS_NTA_C2(2),FUS_NTA_C1(2),Data_type_FUS_NTA);
y_FUS_NTA_050_fit=PN_P_Nstar_fit(x_vec,FUS_NTA_C3(3),FUS_NTA_C2(3),FUS_NTA_C1(3),Data_type_FUS_NTA);
y_FUS_NTA_100_fit=PN_P_Nstar_fit(x_vec,FUS_NTA_C3(4),FUS_NTA_C2(4),FUS_NTA_C1(4),Data_type_FUS_NTA);
y_FUS_NTA_200_fit=PN_P_Nstar_fit(x_vec,FUS_NTA_C3(5),FUS_NTA_C2(5),FUS_NTA_C1(5),Data_type_FUS_NTA);

FUS_NTA_Tension=(FUS_NTA_C2+FUS_NTA_Cm1.*FUS_NTA_C3)/4/pi./R_star_FUS_NTA.^2*KB*300;
FUS_NTA_J0kappa=-FUS_NTA_C1/8/pi./R_star_FUS_NTA*KB*300;
FUS_NTA_shell_core=FUS_NTA_Cm1/3;
FUS_NTA_core_energy_density=FUS_NTA_C3./(4*pi/3*R_star_FUS_NTA.^3);

end


if FUS_DLS
data_DLS_FUS=load('FUS DLS.txt');
Data_type_FUS_DLS='Number';

lower_cluster_radius=20;

Radius_FUS_DLS=data_DLS_FUS(:,1)/2*10^-9;
FUS_DLS_025=data_DLS_FUS(:,2)./sum(data_DLS_FUS(:,2)); % 0.25 micro molar
FUS_DLS_050=data_DLS_FUS(:,4)./sum(data_DLS_FUS(:,4)); % 0.50 micro molar
FUS_DLS_070=data_DLS_FUS(:,6)./sum(data_DLS_FUS(:,6)); % 0.70 micro molar
FUS_DLS_100=data_DLS_FUS(:,8)./sum(data_DLS_FUS(:,8)); % 1.00 micro molar
FUS_DLS_200=data_DLS_FUS(:,10)./sum(data_DLS_FUS(:,10)); % 2.00 micro molar
FUS_DLS_300=data_DLS_FUS(:,12)./sum(data_DLS_FUS(:,12)); % 3.00 micro molar

FUS_DLS_025_error=data_DLS_FUS(:,3)./sum(data_DLS_FUS(:,2)); % 0.25 micro molar
FUS_DLS_050_error=data_DLS_FUS(:,5)./sum(data_DLS_FUS(:,4)); % 0.50 micro molar
FUS_DLS_070_error=data_DLS_FUS(:,7)./sum(data_DLS_FUS(:,6)); % 0.70 micro molar
FUS_DLS_100_error=data_DLS_FUS(:,9)./sum(data_DLS_FUS(:,8)); % 1.00 micro molar
FUS_DLS_200_error=data_DLS_FUS(:,11)./sum(data_DLS_FUS(:,10)); % 2.00 micro molar
FUS_DLS_300_error=data_DLS_FUS(:,13)./sum(data_DLS_FUS(:,12)); % 3.00 micro molar


[~,lower_bound_index]=min(abs(Radius_FUS_DLS-lower_cluster_radius*10^-9));
last_index=length(Radius_FUS_DLS);
% find radius and distribution width
[~,index_FUS_DLS(1)]=max(FUS_DLS_025(lower_bound_index:last_index)); 
[~,index_FUS_DLS(2)]=max(FUS_DLS_050(lower_bound_index:last_index)); 
[~,index_FUS_DLS(3)]=max(FUS_DLS_070(lower_bound_index:last_index));
[~,index_FUS_DLS(4)]=max(FUS_DLS_100(lower_bound_index:last_index)); 
[~,index_FUS_DLS(5)]=max(FUS_DLS_200(lower_bound_index:last_index));
[~,index_FUS_DLS(6)]=max(FUS_DLS_300(lower_bound_index:last_index));

index_FUS_DLS=index_FUS_DLS+lower_bound_index-1;

[R_star_FUS_DLS(1),FUS_DLS_max_vec(1)] =Find_R_Star(Radius_FUS_DLS,FUS_DLS_025,index_FUS_DLS(1));
[R_star_FUS_DLS(2),FUS_DLS_max_vec(2)] =Find_R_Star(Radius_FUS_DLS,FUS_DLS_050,index_FUS_DLS(2));
[R_star_FUS_DLS(3),FUS_DLS_max_vec(3)] =Find_R_Star(Radius_FUS_DLS,FUS_DLS_070,index_FUS_DLS(3));
[R_star_FUS_DLS(4),FUS_DLS_max_vec(4)] =Find_R_Star(Radius_FUS_DLS,FUS_DLS_100,index_FUS_DLS(4));
[R_star_FUS_DLS(5),FUS_DLS_max_vec(5)] =Find_R_Star(Radius_FUS_DLS,FUS_DLS_200,index_FUS_DLS(5));
[R_star_FUS_DLS(6),FUS_DLS_max_vec(6)] =Find_R_Star(Radius_FUS_DLS,FUS_DLS_300,index_FUS_DLS(6));



x_FUS_DLS_025=Radius_FUS_DLS/R_star_FUS_DLS(1);
x_FUS_DLS_050=Radius_FUS_DLS/R_star_FUS_DLS(2);
x_FUS_DLS_070=Radius_FUS_DLS/R_star_FUS_DLS(3);
x_FUS_DLS_100=Radius_FUS_DLS/R_star_FUS_DLS(4);
x_FUS_DLS_200=Radius_FUS_DLS/R_star_FUS_DLS(5);
x_FUS_DLS_300=Radius_FUS_DLS/R_star_FUS_DLS(6);


min_index_FUS_DLS(1)=min_index_finder(x_FUS_DLS_025,index_FUS_DLS(1),x_min);
min_index_FUS_DLS(2)=min_index_finder(x_FUS_DLS_050,index_FUS_DLS(2),x_min);
min_index_FUS_DLS(3)=min_index_finder(x_FUS_DLS_070,index_FUS_DLS(3),x_min);
min_index_FUS_DLS(4)=min_index_finder(x_FUS_DLS_100,index_FUS_DLS(4),x_min);
min_index_FUS_DLS(5)=min_index_finder(x_FUS_DLS_200,index_FUS_DLS(5),x_min);
min_index_FUS_DLS(6)=min_index_finder(x_FUS_DLS_300,index_FUS_DLS(6),x_min);

y_FUS_DLS_025=FUS_DLS_025./FUS_DLS_max_vec(1);
y_FUS_DLS_050=FUS_DLS_050./FUS_DLS_max_vec(2);
y_FUS_DLS_070=FUS_DLS_070./FUS_DLS_max_vec(3);
y_FUS_DLS_100=FUS_DLS_100./FUS_DLS_max_vec(4);
y_FUS_DLS_200=FUS_DLS_200./FUS_DLS_max_vec(5);
y_FUS_DLS_300=FUS_DLS_300./FUS_DLS_max_vec(6);


[~,FUS_DLS_index_last(1)] = Find_Nmax(x_max,Radius_FUS_DLS,y_FUS_DLS_025,y_min,index_FUS_DLS(1));
[~,FUS_DLS_index_last(2)] = Find_Nmax(x_max,Radius_FUS_DLS,y_FUS_DLS_050,y_min,index_FUS_DLS(2));
[~,FUS_DLS_index_last(3)] = Find_Nmax(x_max,Radius_FUS_DLS,y_FUS_DLS_070,y_min,index_FUS_DLS(3));
[~,FUS_DLS_index_last(4)] = Find_Nmax(x_max,Radius_FUS_DLS,y_FUS_DLS_100,y_min,index_FUS_DLS(4));
[~,FUS_DLS_index_last(5)] = Find_Nmax(x_max,Radius_FUS_DLS,y_FUS_DLS_200,y_min,index_FUS_DLS(5));
[~,FUS_DLS_index_last(6)] = Find_Nmax(x_max,Radius_FUS_DLS,y_FUS_DLS_300,y_min,index_FUS_DLS(6));



[FUS_DLS_C3(1),FUS_DLS_C2(1),FUS_DLS_C1(1),FUS_DLS_Cm1(1),FUS_DLS_chi2(1),FUS_DLS_Error_RMS(1)] = LSM_C1C2C3_no_constreint(x_FUS_DLS_025,y_FUS_DLS_025,min_index_FUS_DLS(1),FUS_DLS_index_last(1),6,dC,'log',Data_type_FUS_DLS,0);
[FUS_DLS_C3(2),FUS_DLS_C2(2),FUS_DLS_C1(2),FUS_DLS_Cm1(2),FUS_DLS_chi2(2),FUS_DLS_Error_RMS(2)] = LSM_C1C2C3_no_constreint(x_FUS_DLS_050,y_FUS_DLS_050,min_index_FUS_DLS(2),FUS_DLS_index_last(2),6,dC,'log',Data_type_FUS_DLS,0);
[FUS_DLS_C3(3),FUS_DLS_C2(3),FUS_DLS_C1(3),FUS_DLS_Cm1(3),FUS_DLS_chi2(3),FUS_DLS_Error_RMS(3)] = LSM_C1C2C3_no_constreint(x_FUS_DLS_070,y_FUS_DLS_070,min_index_FUS_DLS(3),FUS_DLS_index_last(3),6,dC,'log',Data_type_FUS_DLS,0);
[FUS_DLS_C3(4),FUS_DLS_C2(4),FUS_DLS_C1(4),FUS_DLS_Cm1(4),FUS_DLS_chi2(4),FUS_DLS_Error_RMS(4)] = LSM_C1C2C3_no_constreint(x_FUS_DLS_100,y_FUS_DLS_100,min_index_FUS_DLS(4),FUS_DLS_index_last(4),6,dC,'log',Data_type_FUS_DLS,0);
[FUS_DLS_C3(5),FUS_DLS_C2(5),FUS_DLS_C1(5),FUS_DLS_Cm1(5),FUS_DLS_chi2(5),FUS_DLS_Error_RMS(5)] = LSM_C1C2C3_no_constreint(x_FUS_DLS_200,y_FUS_DLS_200,min_index_FUS_DLS(5),FUS_DLS_index_last(5),6,dC,'log',Data_type_FUS_DLS,0);
[FUS_DLS_C3(6),FUS_DLS_C2(6),FUS_DLS_C1(6),FUS_DLS_Cm1(6),FUS_DLS_chi2(6),FUS_DLS_Error_RMS(6)] = LSM_C1C2C3_no_constreint(x_FUS_DLS_300,y_FUS_DLS_300,min_index_FUS_DLS(6),FUS_DLS_index_last(6),6,dC,'log',Data_type_FUS_DLS,0);



y_FUS_DLS_025_fit=PN_P_Nstar_fit(x_vec,FUS_DLS_C3(1),FUS_DLS_C2(1),FUS_DLS_C1(1),Data_type_FUS_DLS);
y_FUS_DLS_050_fit=PN_P_Nstar_fit(x_vec,FUS_DLS_C3(2),FUS_DLS_C2(2),FUS_DLS_C1(2),Data_type_FUS_DLS);
y_FUS_DLS_070_fit=PN_P_Nstar_fit(x_vec,FUS_DLS_C3(3),FUS_DLS_C2(3),FUS_DLS_C1(3),Data_type_FUS_DLS);
y_FUS_DLS_100_fit=PN_P_Nstar_fit(x_vec,FUS_DLS_C3(4),FUS_DLS_C2(4),FUS_DLS_C1(4),Data_type_FUS_DLS);
y_FUS_DLS_200_fit=PN_P_Nstar_fit(x_vec,FUS_DLS_C3(5),FUS_DLS_C2(5),FUS_DLS_C1(5),Data_type_FUS_DLS);
y_FUS_DLS_300_fit=PN_P_Nstar_fit(x_vec,FUS_DLS_C3(6),FUS_DLS_C2(6),FUS_DLS_C1(6),Data_type_FUS_DLS);

FUS_DLS_Tension=(FUS_DLS_C2+FUS_DLS_Cm1.*FUS_DLS_C3)/4/pi./R_star_FUS_DLS.^2*KB*300;
FUS_DLS_J0kappa=-FUS_DLS_C1/8/pi./R_star_FUS_DLS*KB*300;
FUS_DLS_shell_core=FUS_DLS_Cm1/3;
FUS_DLS_core_energy_density=FUS_DLS_C3./(4*pi/3*R_star_FUS_DLS.^3);





end


if CPEB4

data_DLS_CPEB4=importdata('CPEB4 DLS.txt');
Data_type_CPEB4='Number';

lower_cluster_radius=1;


Radius_CPEB4=data_DLS_CPEB4(:,1)/2*10^-9;
CPEB4_10=data_DLS_CPEB4(:,2)./sum(data_DLS_CPEB4(:,2)); % 10 micro molar
CPEB4_50=data_DLS_CPEB4(:,4)./sum(data_DLS_CPEB4(:,4)); % 50 micro molar
CPEB4_100=data_DLS_CPEB4(:,6)./sum(data_DLS_CPEB4(:,6)); % 100


[~,lower_bound_index]=min(abs(Radius_CPEB4-lower_cluster_radius*10^-9));
last_index=length(Radius_CPEB4);
% find radius and distribution width
[~,index_CPEB4(1)]=max(CPEB4_10(lower_bound_index:last_index)); 
[~,index_CPEB4(2)]=max(CPEB4_50(lower_bound_index:last_index)); 
[~,index_CPEB4(3)]=max(CPEB4_100(lower_bound_index:last_index));
index_CPEB4=index_CPEB4+lower_bound_index-1;

[R_star_CPEB4(1),CPEB4_max_vec(1)] =Find_R_Star(Radius_CPEB4,CPEB4_10,index_CPEB4(1));
[R_star_CPEB4(2),CPEB4_max_vec(2)] =Find_R_Star(Radius_CPEB4,CPEB4_50,index_CPEB4(2));
[R_star_CPEB4(3),CPEB4_max_vec(3)] =Find_R_Star(Radius_CPEB4,CPEB4_100,index_CPEB4(3));




x_CPEB4_10=Radius_CPEB4/R_star_CPEB4(1);
x_CPEB4_50=Radius_CPEB4/R_star_CPEB4(2);
x_CPEB4_100=Radius_CPEB4/R_star_CPEB4(3);

min_index_CPEB4(1)=min_index_finder(x_CPEB4_10,index_CPEB4(1),x_min);
min_index_CPEB4(2)=min_index_finder(x_CPEB4_50,index_CPEB4(2),x_min);
min_index_CPEB4(3)=min_index_finder(x_CPEB4_100,index_CPEB4(3),x_min);

y_CPEB4_10=CPEB4_10./CPEB4_max_vec(1);
y_CPEB4_50=CPEB4_50./CPEB4_max_vec(2);
y_CPEB4_100=CPEB4_100./CPEB4_max_vec(3);

[~,CPEB4_index_last(1)] = Find_Nmax(x_max,Radius_CPEB4,y_CPEB4_10,y_min,index_CPEB4(1));
[~,CPEB4_index_last(2)] = Find_Nmax(x_max,Radius_CPEB4,y_CPEB4_50,y_min,index_CPEB4(2));
[~,CPEB4_index_last(3)] = Find_Nmax(x_max,Radius_CPEB4,y_CPEB4_100,y_min,index_CPEB4(3));

[CPEB4_C3(1),CPEB4_C2(1),CPEB4_C1(1),CPEB4_Cm1(1),CPEB4_chi2(1),CPEB4_Error_RMS(1)] = LSM_C1C2C3_no_constreint(x_CPEB4_10,y_CPEB4_10,min_index_CPEB4(1),CPEB4_index_last(1),2.5,dC,'log',Data_type_CPEB4,0);
[CPEB4_C3(2),CPEB4_C2(2),CPEB4_C1(2),CPEB4_Cm1(2),CPEB4_chi2(2),CPEB4_Error_RMS(2)] = LSM_C1C2C3_no_constreint(x_CPEB4_50,y_CPEB4_50,min_index_CPEB4(2),CPEB4_index_last(2),2.5,dC,'log',Data_type_CPEB4,0);
[CPEB4_C3(3),CPEB4_C2(3),CPEB4_C1(3),CPEB4_Cm1(3),CPEB4_chi2(3),CPEB4_Error_RMS(3)] = LSM_C1C2C3_no_constreint(x_CPEB4_100,y_CPEB4_100,min_index_CPEB4(3),CPEB4_index_last(3),2.5,dC,'log',Data_type_CPEB4,0);


y_CPEB4_10_fit=PN_P_Nstar_fit(x_vec,CPEB4_C3(1),CPEB4_C2(1),CPEB4_C1(1),Data_type_CPEB4);
y_CPEB4_50_fit=PN_P_Nstar_fit(x_vec,CPEB4_C3(2),CPEB4_C2(2),CPEB4_C1(2),Data_type_CPEB4);
y_CPEB4_100_fit=PN_P_Nstar_fit(x_vec,CPEB4_C3(3),CPEB4_C2(3),CPEB4_C1(3),Data_type_CPEB4);

CPEB4_Tension=(CPEB4_C2+CPEB4_Cm1.*CPEB4_C3)/4/pi./R_star_CPEB4.^2*KB*300;
CPEB4_J0kappa=-CPEB4_C1/8/pi./R_star_CPEB4*KB*300;
CPEB4_shell_core=CPEB4_Cm1/3;
CPEB4_core_energy_density=CPEB4_C3./(4*pi/3*R_star_CPEB4.^3);


end





if hnRNPA3
data_DLS_hnRNPA3=importdata('hnRNPA3.txt');

lower_cluster_radius=1;
Radius_hnRNPA3=data_DLS_hnRNPA3(:,1)/2*10^-9;


Data_type_hnRNPA3='Intensity';
hnRNPA3_0125=data_DLS_hnRNPA3(:,2)./sum(data_DLS_hnRNPA3(:,2));
hnRNPA3_025=data_DLS_hnRNPA3(:,3)./sum(data_DLS_hnRNPA3(:,3));
hnRNPA3_050=data_DLS_hnRNPA3(:,4)./sum(data_DLS_hnRNPA3(:,4));
hnRNPA3_100=data_DLS_hnRNPA3(:,5)./sum(data_DLS_hnRNPA3(:,5)); 
hnRNPA3_200=data_DLS_hnRNPA3(:,6)./sum(data_DLS_hnRNPA3(:,6)); 
hnRNPA3_600=data_DLS_hnRNPA3(:,7)./sum(data_DLS_hnRNPA3(:,7)); 


hnRNPA3_0125=data_DLS_hnRNPA3(:,2)./sum(data_DLS_hnRNPA3(:,2));
hnRNPA3_025=data_DLS_hnRNPA3(:,3)./sum(data_DLS_hnRNPA3(:,3));
hnRNPA3_050=data_DLS_hnRNPA3(:,4)./sum(data_DLS_hnRNPA3(:,4));
hnRNPA3_100=data_DLS_hnRNPA3(:,5)./sum(data_DLS_hnRNPA3(:,5)); 
hnRNPA3_200=data_DLS_hnRNPA3(:,6)./sum(data_DLS_hnRNPA3(:,6)); 
hnRNPA3_600=data_DLS_hnRNPA3(:,7)./sum(data_DLS_hnRNPA3(:,7)); 


[~,lower_bound_index]=min(abs(Radius_hnRNPA3-lower_cluster_radius*10^-9));
last_index=length(Radius_hnRNPA3);
% find radius and distribution width
[~,index_hnRNPA3(1)]=max(hnRNPA3_0125(lower_bound_index:last_index)); 
[~,index_hnRNPA3(2)]=max(hnRNPA3_025(lower_bound_index:last_index)); 
[~,index_hnRNPA3(3)]=max(hnRNPA3_050(lower_bound_index:last_index));
[~,index_hnRNPA3(4)]=max(hnRNPA3_100(lower_bound_index:last_index)); 
[~,index_hnRNPA3(5)]=max(hnRNPA3_200(lower_bound_index:last_index));
[~,index_hnRNPA3(6)]=max(hnRNPA3_600(lower_bound_index:last_index));

index_hnRNPA3=index_hnRNPA3+lower_bound_index-1;

[R_star_hnRNPA3(1),hnRNPA3_max_vec(1)] =Find_R_Star(Radius_hnRNPA3,hnRNPA3_0125,index_hnRNPA3(1));
[R_star_hnRNPA3(2),hnRNPA3_max_vec(2)] =Find_R_Star(Radius_hnRNPA3,hnRNPA3_025,index_hnRNPA3(2));
[R_star_hnRNPA3(3),hnRNPA3_max_vec(3)] =Find_R_Star(Radius_hnRNPA3,hnRNPA3_050,index_hnRNPA3(3));
[R_star_hnRNPA3(4),hnRNPA3_max_vec(4)] =Find_R_Star(Radius_hnRNPA3,hnRNPA3_100,index_hnRNPA3(4));
[R_star_hnRNPA3(5),hnRNPA3_max_vec(5)] =Find_R_Star(Radius_hnRNPA3,hnRNPA3_200,index_hnRNPA3(5));
[R_star_hnRNPA3(6),hnRNPA3_max_vec(6)] =Find_R_Star(Radius_hnRNPA3,hnRNPA3_600,index_hnRNPA3(6));



x_hnRNPA3_0125=Radius_hnRNPA3/R_star_hnRNPA3(1);
x_hnRNPA3_025=Radius_hnRNPA3/R_star_hnRNPA3(2);
x_hnRNPA3_050=Radius_hnRNPA3/R_star_hnRNPA3(3);
x_hnRNPA3_100=Radius_hnRNPA3/R_star_hnRNPA3(4);
x_hnRNPA3_200=Radius_hnRNPA3/R_star_hnRNPA3(5);
x_hnRNPA3_600=Radius_hnRNPA3/R_star_hnRNPA3(6);



min_index_hnRNPA3(1)=min_index_finder(x_hnRNPA3_0125,index_hnRNPA3(1),x_min);
min_index_hnRNPA3(2)=min_index_finder(x_hnRNPA3_025,index_hnRNPA3(2),x_min);
min_index_hnRNPA3(3)=min_index_finder(x_hnRNPA3_050,index_hnRNPA3(3),x_min);
min_index_hnRNPA3(4)=min_index_finder(x_hnRNPA3_100,index_hnRNPA3(4),x_min);
min_index_hnRNPA3(5)=min_index_finder(x_hnRNPA3_200,index_hnRNPA3(5),x_min);
min_index_hnRNPA3(6)=min_index_finder(x_hnRNPA3_600,index_hnRNPA3(6),x_min);


y_hnRNPA3_0125=hnRNPA3_0125./hnRNPA3_max_vec(1);
y_hnRNPA3_025=hnRNPA3_025./hnRNPA3_max_vec(2);
y_hnRNPA3_050=hnRNPA3_050./hnRNPA3_max_vec(3);
y_hnRNPA3_100=hnRNPA3_100./hnRNPA3_max_vec(4);
y_hnRNPA3_200=hnRNPA3_200./hnRNPA3_max_vec(5);
y_hnRNPA3_600=hnRNPA3_600./hnRNPA3_max_vec(6);


[~,hnRNPA3_DLS_index_last(1)] = Find_Nmax(x_max,Radius_hnRNPA3,y_hnRNPA3_0125,y_min,index_hnRNPA3(1));
[~,hnRNPA3_DLS_index_last(2)] = Find_Nmax(x_max,Radius_hnRNPA3,y_hnRNPA3_025,y_min,index_hnRNPA3(2));
[~,hnRNPA3_DLS_index_last(3)] = Find_Nmax(x_max,Radius_hnRNPA3,y_hnRNPA3_050,y_min,index_hnRNPA3(3));
[~,hnRNPA3_DLS_index_last(4)] = Find_Nmax(x_max,Radius_hnRNPA3,y_hnRNPA3_100,y_min,index_hnRNPA3(4));
[~,hnRNPA3_DLS_index_last(5)] = Find_Nmax(x_max,Radius_hnRNPA3,y_hnRNPA3_200,y_min,index_hnRNPA3(5));
[~,hnRNPA3_DLS_index_last(6)] = Find_Nmax(x_max,Radius_hnRNPA3,y_hnRNPA3_600,y_min,index_hnRNPA3(6));


[hnRNPA3_C3(1),hnRNPA3_C2(1),hnRNPA3_C1(1),hnRNPA3_DLS_Cm1(1),hnRNPA3_chi2(1),hnRNPA3_Error_RMS(1)] = LSM_C1C2C3_no_constreint(x_hnRNPA3_0125,y_hnRNPA3_0125,min_index_hnRNPA3(1),hnRNPA3_DLS_index_last(1),6,dC,'log',Data_type_hnRNPA3,0);
[hnRNPA3_C3(2),hnRNPA3_C2(2),hnRNPA3_C1(2),hnRNPA3_DLS_Cm1(2),hnRNPA3_chi2(2),hnRNPA3_Error_RMS(2)] = LSM_C1C2C3_no_constreint(x_hnRNPA3_025,y_hnRNPA3_025,min_index_hnRNPA3(2),hnRNPA3_DLS_index_last(2),6,dC,'log',Data_type_hnRNPA3,0);
[hnRNPA3_C3(3),hnRNPA3_C2(3),hnRNPA3_C1(3),hnRNPA3_DLS_Cm1(3),hnRNPA3_chi2(3),hnRNPA3_Error_RMS(3)] = LSM_C1C2C3_no_constreint(x_hnRNPA3_050,y_hnRNPA3_050,min_index_hnRNPA3(3),hnRNPA3_DLS_index_last(3),6,dC,'log',Data_type_hnRNPA3,0);
[hnRNPA3_C3(4),hnRNPA3_C2(4),hnRNPA3_C1(4),hnRNPA3_DLS_Cm1(4),hnRNPA3_chi2(4),hnRNPA3_Error_RMS(4)] = LSM_C1C2C3_no_constreint(x_hnRNPA3_100,y_hnRNPA3_100,min_index_hnRNPA3(4),hnRNPA3_DLS_index_last(4),6,dC,'log',Data_type_hnRNPA3,0);
[hnRNPA3_C3(5),hnRNPA3_C2(5),hnRNPA3_C1(5),hnRNPA3_DLS_Cm1(5),hnRNPA3_chi2(5),hnRNPA3_Error_RMS(5)] = LSM_C1C2C3_no_constreint(x_hnRNPA3_200,y_hnRNPA3_200,min_index_hnRNPA3(5),hnRNPA3_DLS_index_last(5),6,dC,'log',Data_type_hnRNPA3,0);
[hnRNPA3_C3(6),hnRNPA3_C2(6),hnRNPA3_C1(6),hnRNPA3_DLS_Cm1(6),hnRNPA3_chi2(6),hnRNPA3_Error_RMS(6)] = LSM_C1C2C3_no_constreint(x_hnRNPA3_600,y_hnRNPA3_600,min_index_hnRNPA3(6),hnRNPA3_DLS_index_last(6),6,dC,'log',Data_type_hnRNPA3,0);



y_hnRNPA3_0125_fit=PN_P_Nstar_fit(x_vec,hnRNPA3_C3(1),hnRNPA3_C2(1),hnRNPA3_C1(1),Data_type_hnRNPA3);
y_hnRNPA3_025_fit=PN_P_Nstar_fit(x_vec,hnRNPA3_C3(2),hnRNPA3_C2(2),hnRNPA3_C1(2),Data_type_hnRNPA3);
y_hnRNPA3_050_fit=PN_P_Nstar_fit(x_vec,hnRNPA3_C3(3),hnRNPA3_C2(3),hnRNPA3_C1(3),Data_type_hnRNPA3);
y_hnRNPA3_100_fit=PN_P_Nstar_fit(x_vec,hnRNPA3_C3(4),hnRNPA3_C2(4),hnRNPA3_C1(4),Data_type_hnRNPA3);
y_hnRNPA3_200_fit=PN_P_Nstar_fit(x_vec,hnRNPA3_C3(5),hnRNPA3_C2(5),hnRNPA3_C1(5),Data_type_hnRNPA3);
y_hnRNPA3_600_fit=PN_P_Nstar_fit(x_vec,hnRNPA3_C3(6),hnRNPA3_C2(6),hnRNPA3_C1(6),Data_type_hnRNPA3);

hnRNPA3_DLS_Tension=(hnRNPA3_C2+hnRNPA3_DLS_Cm1.*hnRNPA3_C3)/4/pi./R_star_hnRNPA3.^2*KB*300;
hnRNPA3_DLS_J0kappa=-hnRNPA3_C1/8/pi./R_star_hnRNPA3*KB*300;
hnRNPA3_DLS_shell_core=hnRNPA3_DLS_Cm1/3;
hnRNPA3_DLS_core_energy_density=hnRNPA3_C3./(4*pi/3*R_star_hnRNPA3.^3);


y_hnRNPA3_0125_fit=y_hnRNPA3_0125_fit./x_vec.^6;
y_hnRNPA3_025_fit=y_hnRNPA3_025_fit./x_vec.^6;
y_hnRNPA3_050_fit=y_hnRNPA3_050_fit./x_vec.^6;
y_hnRNPA3_100_fit=y_hnRNPA3_100_fit./x_vec.^6;
y_hnRNPA3_200_fit=y_hnRNPA3_200_fit./x_vec.^6;
y_hnRNPA3_600_fit=y_hnRNPA3_600_fit./x_vec.^6;

y_hnRNPA3_0125=y_hnRNPA3_0125./x_hnRNPA3_0125.^6;
y_hnRNPA3_025=y_hnRNPA3_025./x_hnRNPA3_025.^6;
y_hnRNPA3_050=y_hnRNPA3_050./x_hnRNPA3_050.^6;
y_hnRNPA3_100=y_hnRNPA3_100./x_hnRNPA3_100.^6;
y_hnRNPA3_200=y_hnRNPA3_200./x_hnRNPA3_200.^6;
y_hnRNPA3_600=y_hnRNPA3_600./x_hnRNPA3_600.^6;


end


if EWSR1
data_DLS_EWSR1=importdata("EWSR1.txt");


Data_type_EWSR1='Intensity';

lower_cluster_radius=1;

Radius_EWSR1=data_DLS_EWSR1(:,1)/2*10^-9;
EWSR1_025=data_DLS_EWSR1(:,2)./sum(data_DLS_EWSR1(:,2));
EWSR1_050=data_DLS_EWSR1(:,3)./sum(data_DLS_EWSR1(:,3));
EWSR1_100=data_DLS_EWSR1(:,4)./sum(data_DLS_EWSR1(:,4)); 
EWSR1_200=data_DLS_EWSR1(:,5)./sum(data_DLS_EWSR1(:,5)); 


[~,lower_bound_index]=min(abs(Radius_EWSR1-lower_cluster_radius*10^-9));
last_index=length(Radius_EWSR1);
% find radius and distribution width
[~,index_EWSR1(1)]=max(EWSR1_025(lower_bound_index:last_index)); 
[~,index_EWSR1(2)]=max(EWSR1_050(lower_bound_index:last_index));
[~,index_EWSR1(3)]=max(EWSR1_100(lower_bound_index:last_index)); 
[~,index_EWSR1(4)]=max(EWSR1_200(lower_bound_index:last_index));

index_EWSR1=index_EWSR1+lower_bound_index-1;

[R_star_EWSR1(1),EWSR1_max_vec(1)] =Find_R_Star(Radius_EWSR1,EWSR1_025,index_EWSR1(1));
[R_star_EWSR1(2),EWSR1_max_vec(2)] =Find_R_Star(Radius_EWSR1,EWSR1_050,index_EWSR1(2));
[R_star_EWSR1(3),EWSR1_max_vec(3)] =Find_R_Star(Radius_EWSR1,EWSR1_100,index_EWSR1(3));
[R_star_EWSR1(4),EWSR1_max_vec(4)] =Find_R_Star(Radius_EWSR1,EWSR1_200,index_EWSR1(4));



x_EWSR1_025=Radius_EWSR1/R_star_EWSR1(1);
x_EWSR1_050=Radius_EWSR1/R_star_EWSR1(2);
x_EWSR1_100=Radius_EWSR1/R_star_EWSR1(3);
x_EWSR1_200=Radius_EWSR1/R_star_EWSR1(4);



min_index_EWSR1(1)=min_index_finder(x_EWSR1_025,index_EWSR1(1),x_min);
min_index_EWSR1(2)=min_index_finder(x_EWSR1_050,index_EWSR1(2),x_min);
min_index_EWSR1(3)=min_index_finder(x_EWSR1_100,index_EWSR1(3),x_min);
min_index_EWSR1(4)=min_index_finder(x_EWSR1_200,index_EWSR1(4),x_min);


y_EWSR1_025=EWSR1_025./EWSR1_max_vec(1);
y_EWSR1_050=EWSR1_050./EWSR1_max_vec(2);
y_EWSR1_100=EWSR1_100./EWSR1_max_vec(3);
y_EWSR1_200=EWSR1_200./EWSR1_max_vec(4);


[~,EWSR1_index_last(1)] = Find_Nmax(x_max,Radius_EWSR1,y_EWSR1_025,y_min,index_EWSR1(1));
[~,EWSR1_index_last(2)] = Find_Nmax(x_max,Radius_EWSR1,y_EWSR1_050,y_min,index_EWSR1(2));
[~,EWSR1_index_last(3)] = Find_Nmax(x_max,Radius_EWSR1,y_EWSR1_100,y_min,index_EWSR1(3));
[~,EWSR1_index_last(4)] = Find_Nmax(x_max,Radius_EWSR1,y_EWSR1_200,y_min,index_EWSR1(4));


[EWSR1_C3(1),EWSR1_C2(1),EWSR1_C1(1),EWSR1_DLS_Cm1(1),EWSR1_chi2(1),EWSR1_Error_RMS(1)] = LSM_C1C2C3_no_constreint(x_EWSR1_025,y_EWSR1_025,min_index_EWSR1(1),EWSR1_index_last(1),6,dC,'log',Data_type_EWSR1,0);
[EWSR1_C3(2),EWSR1_C2(2),EWSR1_C1(2),EWSR1_DLS_Cm1(2),EWSR1_chi2(2),EWSR1_Error_RMS(2)] = LSM_C1C2C3_no_constreint(x_EWSR1_050,y_EWSR1_050,min_index_EWSR1(2),EWSR1_index_last(2),6,dC,'log',Data_type_EWSR1,0);
[EWSR1_C3(3),EWSR1_C2(3),EWSR1_C1(3),EWSR1_DLS_Cm1(3),EWSR1_chi2(3),EWSR1_Error_RMS(3)] = LSM_C1C2C3_no_constreint(x_EWSR1_100,y_EWSR1_100,min_index_EWSR1(3),EWSR1_index_last(3),6,dC,'log',Data_type_EWSR1,0);
[EWSR1_C3(4),EWSR1_C2(4),EWSR1_C1(4),EWSR1_DLS_Cm1(4),EWSR1_chi2(4),EWSR1_Error_RMS(4)] = LSM_C1C2C3_no_constreint(x_EWSR1_200,y_EWSR1_200,min_index_EWSR1(4),EWSR1_index_last(4),6,dC,'log',Data_type_EWSR1,0);

y_EWSR1_025_fit=PN_P_Nstar_fit(x_vec,EWSR1_C3(1),EWSR1_C2(1),EWSR1_C1(1),Data_type_EWSR1);
y_EWSR1_050_fit=PN_P_Nstar_fit(x_vec,EWSR1_C3(2),EWSR1_C2(2),EWSR1_C1(2),Data_type_EWSR1);
y_EWSR1_100_fit=PN_P_Nstar_fit(x_vec,EWSR1_C3(3),EWSR1_C2(3),EWSR1_C1(3),Data_type_EWSR1);
y_EWSR1_200_fit=PN_P_Nstar_fit(x_vec,EWSR1_C3(4),EWSR1_C2(4),EWSR1_C1(4),Data_type_EWSR1);

EWSR1_DLS_Tension=(EWSR1_C2+EWSR1_DLS_Cm1.*EWSR1_C3)/4/pi./R_star_EWSR1.^2*KB*300;
EWSR1_DLS_J0kappa=-EWSR1_C1/8/pi./R_star_EWSR1*KB*300;
EWSR1_DLS_shell_core=EWSR1_DLS_Cm1/3;
EWSR1_DLS_core_energy_density=EWSR1_C3./(4*pi/3*R_star_EWSR1.^3);


y_EWSR1_025_fit=y_EWSR1_025_fit./x_vec.^6;
y_EWSR1_050_fit=y_EWSR1_050_fit./x_vec.^6;
y_EWSR1_100_fit=y_EWSR1_100_fit./x_vec.^6;
y_EWSR1_200_fit=y_EWSR1_200_fit./x_vec.^6;

y_EWSR1_025=y_EWSR1_025./x_EWSR1_025.^6;
y_EWSR1_050=y_EWSR1_050./x_EWSR1_050.^6;
y_EWSR1_100=y_EWSR1_100./x_EWSR1_100.^6;
y_EWSR1_200=y_EWSR1_200./x_EWSR1_200.^6;


end

if TAF15 
data_DLS_TAF15=importdata('TAF15.txt');

Data_type_TAF15='Intensity';

lower_cluster_radius=1;

Radius_TAF15=data_DLS_TAF15(:,1)/2*10^-9;
TAF15_0125=data_DLS_TAF15(:,2)./sum(data_DLS_TAF15(:,2));
TAF15_025=data_DLS_TAF15(:,3)./sum(data_DLS_TAF15(:,3));
TAF15_050=data_DLS_TAF15(:,4)./sum(data_DLS_TAF15(:,4));
TAF15_100=data_DLS_TAF15(:,5)./sum(data_DLS_TAF15(:,5)); 
TAF15_200=data_DLS_TAF15(:,6)./sum(data_DLS_TAF15(:,6)); 


[~,lower_bound_index]=min(abs(Radius_TAF15-lower_cluster_radius*10^-9));
last_index=length(Radius_TAF15);
% find radius and distribution width
[~,index_TAF15(1)]=max(TAF15_0125(lower_bound_index:last_index)); 
[~,index_TAF15(2)]=max(TAF15_025(lower_bound_index:last_index)); 
[~,index_TAF15(3)]=max(TAF15_050(lower_bound_index:last_index));
[~,index_TAF15(4)]=max(TAF15_100(lower_bound_index:last_index)); 
[~,index_TAF15(5)]=max(TAF15_200(lower_bound_index:last_index));

index_TAF15=index_TAF15+lower_bound_index-1;

[R_star_TAF15(1),TAF15_max_vec(1)] =Find_R_Star(Radius_TAF15,TAF15_0125,index_TAF15(1));
[R_star_TAF15(2),TAF15_max_vec(2)] =Find_R_Star(Radius_TAF15,TAF15_025,index_TAF15(2));
[R_star_TAF15(3),TAF15_max_vec(3)] =Find_R_Star(Radius_TAF15,TAF15_050,index_TAF15(3));
[R_star_TAF15(4),TAF15_max_vec(4)] =Find_R_Star(Radius_TAF15,TAF15_100,index_TAF15(4));
[R_star_TAF15(5),TAF15_max_vec(5)] =Find_R_Star(Radius_TAF15,TAF15_200,index_TAF15(5));



x_TAF15_0125=Radius_TAF15/R_star_TAF15(1);
x_TAF15_025=Radius_TAF15/R_star_TAF15(2);
x_TAF15_050=Radius_TAF15/R_star_TAF15(3);
x_TAF15_100=Radius_TAF15/R_star_TAF15(4);
x_TAF15_200=Radius_TAF15/R_star_TAF15(5);



min_index_TAF15(1)=min_index_finder(x_TAF15_0125,index_TAF15(1),x_min);
min_index_TAF15(2)=min_index_finder(x_TAF15_025,index_TAF15(2),x_min);
min_index_TAF15(3)=min_index_finder(x_TAF15_050,index_TAF15(3),x_min);
min_index_TAF15(4)=min_index_finder(x_TAF15_100,index_TAF15(4),x_min);
min_index_TAF15(5)=min_index_finder(x_TAF15_200,index_TAF15(5),x_min);


y_TAF15_0125=TAF15_0125./TAF15_max_vec(1);
y_TAF15_025=TAF15_025./TAF15_max_vec(2);
y_TAF15_050=TAF15_050./TAF15_max_vec(3);
y_TAF15_100=TAF15_100./TAF15_max_vec(4);
y_TAF15_200=TAF15_200./TAF15_max_vec(5);


[~,TAF15_DLS_index_last(1)] = Find_Nmax(x_max,Radius_TAF15,y_TAF15_0125,y_min,index_TAF15(1));
[~,TAF15_DLS_index_last(2)] = Find_Nmax(x_max,Radius_TAF15,y_TAF15_025,y_min,index_TAF15(2));
[~,TAF15_DLS_index_last(3)] = Find_Nmax(x_max,Radius_TAF15,y_TAF15_050,y_min,index_TAF15(3));
[~,TAF15_DLS_index_last(4)] = Find_Nmax(x_max,Radius_TAF15,y_TAF15_100,y_min,index_TAF15(4));
[~,TAF15_DLS_index_last(5)] = Find_Nmax(x_max,Radius_TAF15,y_TAF15_200,y_min,index_TAF15(5));


[TAF15_C3(1),TAF15_C2(1),TAF15_C1(1),TAF15_DLS_Cm1(1),TAF15_chi2(1),TAF15_Error_RMS(1)] = LSM_C1C2C3_no_constreint(x_TAF15_0125,y_TAF15_0125,min_index_TAF15(1),TAF15_DLS_index_last(1),6,dC,'log',Data_type_TAF15,0);
[TAF15_C3(2),TAF15_C2(2),TAF15_C1(2),TAF15_DLS_Cm1(2),TAF15_chi2(2),TAF15_Error_RMS(2)] = LSM_C1C2C3_no_constreint(x_TAF15_025,y_TAF15_025,min_index_TAF15(2),TAF15_DLS_index_last(2),6,dC,'log',Data_type_TAF15,0);
[TAF15_C3(3),TAF15_C2(3),TAF15_C1(3),TAF15_DLS_Cm1(3),TAF15_chi2(3),TAF15_Error_RMS(3)] = LSM_C1C2C3_no_constreint(x_TAF15_050,y_TAF15_050,min_index_TAF15(3),TAF15_DLS_index_last(3),6,dC,'log',Data_type_TAF15,0);
[TAF15_C3(4),TAF15_C2(4),TAF15_C1(4),TAF15_DLS_Cm1(4),TAF15_chi2(4),TAF15_Error_RMS(4)] = LSM_C1C2C3_no_constreint(x_TAF15_100,y_TAF15_100,min_index_TAF15(4),TAF15_DLS_index_last(4),6,dC,'log',Data_type_TAF15,0);
[TAF15_C3(5),TAF15_C2(5),TAF15_C1(5),TAF15_DLS_Cm1(5),TAF15_chi2(5),TAF15_Error_RMS(5)] = LSM_C1C2C3_no_constreint(x_TAF15_200,y_TAF15_200,min_index_TAF15(5),TAF15_DLS_index_last(5),6,dC,'log',Data_type_TAF15,0);




y_TAF15_0125_fit=PN_P_Nstar_fit(x_vec,TAF15_C3(1),TAF15_C2(1),TAF15_C1(1),Data_type_TAF15);
y_TAF15_025_fit=PN_P_Nstar_fit(x_vec,TAF15_C3(2),TAF15_C2(2),TAF15_C1(2),Data_type_TAF15);
y_TAF15_050_fit=PN_P_Nstar_fit(x_vec,TAF15_C3(3),TAF15_C2(3),TAF15_C1(3),Data_type_TAF15);
y_TAF15_100_fit=PN_P_Nstar_fit(x_vec,TAF15_C3(4),TAF15_C2(4),TAF15_C1(4),Data_type_TAF15);
y_TAF15_200_fit=PN_P_Nstar_fit(x_vec,TAF15_C3(5),TAF15_C2(5),TAF15_C1(5),Data_type_TAF15);


y_TAF15_0125_fit=y_TAF15_0125_fit./x_vec.^6;
y_TAF15_025_fit=y_TAF15_025_fit./x_vec.^6;
y_TAF15_050_fit=y_TAF15_050_fit./x_vec.^6;
y_TAF15_100_fit=y_TAF15_100_fit./x_vec.^6;
y_TAF15_200_fit=y_TAF15_200_fit./x_vec.^6;

y_TAF15_0125=y_TAF15_0125./x_TAF15_0125.^6;
y_TAF15_025=y_TAF15_025./x_TAF15_025.^6;
y_TAF15_050=y_TAF15_050./x_TAF15_050.^6;
y_TAF15_100=y_TAF15_100./x_TAF15_100.^6;
y_TAF15_200=y_TAF15_200./x_TAF15_200.^6;



TAF15_DLS_Tension=(TAF15_C2+TAF15_DLS_Cm1.*TAF15_C3)/4/pi./R_star_TAF15.^2*KB*300;
TAF15_DLS_J0kappa=-TAF15_C1/8/pi./R_star_TAF15*KB*300;
TAF15_DLS_shell_core=TAF15_DLS_Cm1/3;
TAF15_DLS_core_energy_density=TAF15_C3./(4*pi/3*R_star_TAF15.^3);

end
%% plot fit
figure(5)
hold on
x_test=linspace(0,x_max,1000);
%FUS
if FUS_NTA

s1=scatter(x_FUS_NTA_0125,y_FUS_NTA_0125,72,'MarkerEdgeColor',Col_0125,'LineWidth',2); set(get(get(s1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s2=scatter(x_FUS_NTA_025,y_FUS_NTA_025,72,'MarkerEdgeColor',Col_025,'LineWidth',2); set(get(get(s2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s3=scatter(x_FUS_NTA_050,y_FUS_NTA_050,72,'MarkerEdgeColor',Col_050,'LineWidth',2); set(get(get(s3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s4=scatter(x_FUS_NTA_100,y_FUS_NTA_100,72,'MarkerEdgeColor',Col_100,'LineWidth',2); set(get(get(s4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s5=scatter(x_FUS_NTA_200,y_FUS_NTA_200,72,'MarkerEdgeColor',Col_200,'LineWidth',2); set(get(get(s5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

plot(x_vec,y_FUS_NTA_0125_fit,'--','Color',Col_0125,'LineWidth',2,'DisplayName',['0.125 {\mu}M']);
plot(x_vec,y_FUS_NTA_025_fit,'--','Color',Col_025,'LineWidth',2,'DisplayName','0.25 {\mu}M');
plot(x_vec,y_FUS_NTA_050_fit,'--','Color',Col_050,'LineWidth',2,'DisplayName','0.50 {\mu}M');
plot(x_vec,y_FUS_NTA_100_fit,'--','Color',Col_100,'LineWidth',2,'DisplayName',['1.00 {\mu}M']);
plot(x_vec,y_FUS_NTA_200_fit,'--','Color',Col_200,'LineWidth',2,'DisplayName',['2.00 {\mu}M']);

end

if FUS_DLS

s1=scatter(x_FUS_DLS_025,y_FUS_DLS_025,72,'MarkerEdgeColor',Col_025,'LineWidth',2); set(get(get(s1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s2=scatter(x_FUS_DLS_050,y_FUS_DLS_050,72,'MarkerEdgeColor',Col_050,'LineWidth',2); set(get(get(s2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s3=scatter(x_FUS_DLS_070,y_FUS_DLS_070,72,'MarkerEdgeColor',Col_070,'LineWidth',2); set(get(get(s3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s4=scatter(x_FUS_DLS_100,y_FUS_DLS_100,72,'MarkerEdgeColor',Col_100,'LineWidth',2); set(get(get(s4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s5=scatter(x_FUS_DLS_200,y_FUS_DLS_200,72,'MarkerEdgeColor',Col_200,'LineWidth',2); set(get(get(s5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s6=scatter(x_FUS_DLS_300,y_FUS_DLS_300,72,'MarkerEdgeColor',Col_300,'LineWidth',2); set(get(get(s6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


plot(x_vec,y_FUS_DLS_025_fit,'--','Color',Col_025,'LineWidth',2,'DisplayName','0.25 {\mu}M');
plot(x_vec,y_FUS_DLS_050_fit,'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.50 {\mu}M']);
plot(x_vec,y_FUS_DLS_070_fit,'--','Color',Col_070,'LineWidth',2,'DisplayName',['0.70 {\mu}M']);
plot(x_vec,y_FUS_DLS_100_fit,'--','Color',Col_100,'LineWidth',2,'DisplayName',['1.00 {\mu}M']);
plot(x_vec,y_FUS_DLS_200_fit,'--','Color',Col_200,'LineWidth',2,'DisplayName',['2.00 {\mu}M']);
plot(x_vec,y_FUS_DLS_300_fit,'--','Color',Col_300,'LineWidth',2,'DisplayName',['3.00 {\mu}M']);

end
%CPEB4 conc
if CPEB4
s6=scatter(x_CPEB4_10,y_CPEB4_10,72,'MarkerEdgeColor',Col_025,'LineWidth',2);   set(get(get(s6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s7=scatter(x_CPEB4_50,y_CPEB4_50,72,'MarkerEdgeColor',Col_050,'LineWidth',2);   set(get(get(s7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
s8=scatter(x_CPEB4_100,y_CPEB4_100,72,'MarkerEdgeColor',Col_100,'LineWidth',2); set(get(get(s8,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 

plot(x_vec,y_CPEB4_10_fit,'--','Color',Col_025,'LineWidth',2,'DisplayName',['10 {\mu}M']);
plot(x_vec,y_CPEB4_50_fit,'--','Color',Col_050,'LineWidth',2,'DisplayName',['50 {\mu}M']);
plot(x_vec,y_CPEB4_100_fit,'--','Color',Col_100,'LineWidth',2,'DisplayName',['100 {\mu}M']);

end



%hnRNPA3
if hnRNPA3
s15=scatter(x_hnRNPA3_0125,y_hnRNPA3_0125,72,'MarkerEdgeColor',Col_0125,'LineWidth',2);
s16=scatter(x_hnRNPA3_025,y_hnRNPA3_025,72,'MarkerEdgeColor',Col_025,'LineWidth',2);
s17=scatter(x_hnRNPA3_050,y_hnRNPA3_050,72,'MarkerEdgeColor',Col_050,'LineWidth',2);
s18=scatter(x_hnRNPA3_100,y_hnRNPA3_100,72,'MarkerEdgeColor',Col_100,'LineWidth',2);
s19=scatter(x_hnRNPA3_200,y_hnRNPA3_200,72,'MarkerEdgeColor',Col_200,'LineWidth',2);
s20=scatter(x_hnRNPA3_600,y_hnRNPA3_600,72,'MarkerEdgeColor',Col_600,'LineWidth',2);

plot(x_vec,y_hnRNPA3_0125_fit,'--','Color',Col_0125,'LineWidth',2,'DisplayName',['0.125 {\mu}M']);
plot(x_vec,y_hnRNPA3_025_fit,'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,y_hnRNPA3_050_fit,'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.5 {\mu}M']);
plot(x_vec,y_hnRNPA3_100_fit,'--','Color',Col_100,'LineWidth',2,'DisplayName',['1 {\mu}M']);
plot(x_vec,y_hnRNPA3_200_fit,'--','Color',Col_200,'LineWidth',2,'DisplayName',['2 {\mu}M']);
plot(x_vec,y_hnRNPA3_600_fit,'--','Color',Col_600,'LineWidth',2,'DisplayName',['6 {\mu}M']);

set(get(get(s15,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s16,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s17,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s18,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s19,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s20,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

 
end

%EWSR1
if EWSR1

s16=scatter(x_EWSR1_025,y_EWSR1_025,72,'MarkerEdgeColor',Col_025,'LineWidth',2);
s17=scatter(x_EWSR1_050,y_EWSR1_050,72,'MarkerEdgeColor',Col_050,'LineWidth',2);
s18=scatter(x_EWSR1_100,y_EWSR1_100,72,'MarkerEdgeColor',Col_100,'LineWidth',2);
s19=scatter(x_EWSR1_200,y_EWSR1_200,72,'MarkerEdgeColor',Col_200,'LineWidth',2);

plot(x_vec,y_EWSR1_025_fit,'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,y_EWSR1_050_fit,'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.5 {\mu}M']);
plot(x_vec,y_EWSR1_100_fit,'--','Color',Col_100,'LineWidth',2,'DisplayName',['1 {\mu}M']);
plot(x_vec,y_EWSR1_200_fit,'--','Color',Col_200,'LineWidth',2,'DisplayName',['2 {\mu}M']);


set(get(get(s16,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s17,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s18,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s19,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

%TAF15
if TAF15
s15=scatter(x_TAF15_0125,y_TAF15_0125,72,'MarkerEdgeColor',Col_0125,'LineWidth',2);
s16=scatter(x_TAF15_025,y_TAF15_025,72,'MarkerEdgeColor',Col_025,'LineWidth',2);
s17=scatter(x_TAF15_050,y_TAF15_050,72,'MarkerEdgeColor',Col_050,'LineWidth',2);
s18=scatter(x_TAF15_100,y_TAF15_100,72,'MarkerEdgeColor',Col_100,'LineWidth',2);
s19=scatter(x_TAF15_200,y_TAF15_200,72,'MarkerEdgeColor',Col_200,'LineWidth',2);

plot(x_vec,y_TAF15_0125_fit,'--','Color',Col_0125,'LineWidth',2,'DisplayName','0.125 {\mu}M');
plot(x_vec,y_TAF15_025_fit,'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,y_TAF15_050_fit,'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.5 {\mu}M']);
plot(x_vec,y_TAF15_100_fit,'--','Color',Col_100,'LineWidth',2,'DisplayName',['1 {\mu}M']);
plot(x_vec,y_TAF15_200_fit,'--','Color',Col_200,'LineWidth',2,'DisplayName',['2 {\mu}M']);

set(get(get(s15,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s16,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s17,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s18,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s19,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

p1=plot([1 1],[0 1],'--k','LineWidth',1); set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p2=plot([0 y_min],[3 y_min],'--k','LineWidth',1); set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
xL=xlabel('R/R^*');
yL=ylabel('P_R/P_{R^*}');
set(gca,'FontSize',18);
set(gca,'FontWeight','bold');
set(gca,'LineWidth',1);
xL.FontSize=18;
xL.FontWeight='bold';
yL.FontSize=18;
yL.FontWeight='bold';
set(gcf,'color','white');
xlim([x_min x_max]);
ylim([y_min 1])

legend('FontSize',10)


%% plot fit
figure(11)
hold on
x_test=linspace(0,x_max,1000);
%FUS
if FUS_NTA

s1=scatter(x_FUS_NTA_0125,log(y_FUS_NTA_0125),72,'MarkerEdgeColor',Col_0125,'LineWidth',2); set(get(get(s1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s2=scatter(x_FUS_NTA_025,log(y_FUS_NTA_025),72,'MarkerEdgeColor',Col_025,'LineWidth',2); set(get(get(s2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s3=scatter(x_FUS_NTA_050,log(y_FUS_NTA_050),72,'MarkerEdgeColor',Col_050,'LineWidth',2); set(get(get(s3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s4=scatter(x_FUS_NTA_100,log(y_FUS_NTA_100),72,'MarkerEdgeColor',Col_100,'LineWidth',2); set(get(get(s4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s5=scatter(x_FUS_NTA_200,log(y_FUS_NTA_200),72,'MarkerEdgeColor',Col_200,'LineWidth',2); set(get(get(s5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

plot(x_vec,log(y_FUS_NTA_0125_fit),'--','Color',Col_0125,'LineWidth',2,'DisplayName',['0.125 {\mu}M']);
plot(x_vec,log(y_FUS_NTA_025_fit),'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,log(y_FUS_NTA_050_fit),'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.50 {\mu}M']);
plot(x_vec,log(y_FUS_NTA_100_fit),'--','Color',Col_100,'LineWidth',2,'DisplayName',['1.00 {\mu}M']);
plot(x_vec,log(y_FUS_NTA_200_fit),'--','Color',Col_200,'LineWidth',2,'DisplayName',['2.00 {\mu}M']);

end

if FUS_DLS

s1=scatter(x_FUS_DLS_025,log(y_FUS_DLS_025),72,'MarkerEdgeColor',Col_025,'LineWidth',2); set(get(get(s1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s2=scatter(x_FUS_DLS_050,log(y_FUS_DLS_050),72,'MarkerEdgeColor',Col_050,'LineWidth',2); set(get(get(s2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s3=scatter(x_FUS_DLS_070,log(y_FUS_DLS_070),72,'MarkerEdgeColor',Col_070,'LineWidth',2); set(get(get(s3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s4=scatter(x_FUS_DLS_100,log(y_FUS_DLS_100),72,'MarkerEdgeColor',Col_100,'LineWidth',2); set(get(get(s4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s5=scatter(x_FUS_DLS_200,log(y_FUS_DLS_200),72,'MarkerEdgeColor',Col_200,'LineWidth',2); set(get(get(s5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s6=scatter(x_FUS_DLS_300,log(y_FUS_DLS_300),72,'MarkerEdgeColor',Col_300,'LineWidth',2); set(get(get(s6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


plot(x_vec,log(y_FUS_DLS_025_fit),'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,log(y_FUS_DLS_050_fit),'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.50 {\mu}M']);
plot(x_vec,log(y_FUS_DLS_070_fit),'--','Color',Col_070,'LineWidth',2,'DisplayName',['0.70 {\mu}M']);
plot(x_vec,log(y_FUS_DLS_100_fit),'--','Color',Col_100,'LineWidth',2,'DisplayName',['1.00 {\mu}M']);
plot(x_vec,log(y_FUS_DLS_200_fit),'--','Color',Col_200,'LineWidth',2,'DisplayName',['2.00 {\mu}M']);
plot(x_vec,log(y_FUS_DLS_300_fit),'--','Color',Col_300,'LineWidth',2,'DisplayName',['3.00 {\mu}M']);


end
%CPEB4 conc
if CPEB4
s6=scatter(x_CPEB4_10,log(y_CPEB4_10),72,'MarkerEdgeColor',Col_025,'LineWidth',2);   set(get(get(s6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s7=scatter(x_CPEB4_50,log(y_CPEB4_50),72,'MarkerEdgeColor',Col_050,'LineWidth',2);   set(get(get(s7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
s8=scatter(x_CPEB4_100,log(y_CPEB4_100),72,'MarkerEdgeColor',Col_100,'LineWidth',2); set(get(get(s8,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 

plot(x_vec,log(y_CPEB4_10_fit),'--','Color',Col_025,'LineWidth',2,'DisplayName',['10 {\mu}M']);
plot(x_vec,log(y_CPEB4_50_fit),'--','Color',Col_050,'LineWidth',2,'DisplayName',['50 {\mu}M']);
plot(x_vec,log(y_CPEB4_100_fit),'--','Color',Col_100,'LineWidth',2,'DisplayName',['100 {\mu}M']);

end


%hnRNPA3
if hnRNPA3
s15=scatter(x_hnRNPA3_0125,log(y_hnRNPA3_0125),72,'MarkerEdgeColor',Col_0125,'LineWidth',2);
s16=scatter(x_hnRNPA3_025,log(y_hnRNPA3_025),72,'MarkerEdgeColor',Col_025,'LineWidth',2);
s17=scatter(x_hnRNPA3_050,log(y_hnRNPA3_050),72,'MarkerEdgeColor',Col_050,'LineWidth',2);
s18=scatter(x_hnRNPA3_100,log(y_hnRNPA3_100),72,'MarkerEdgeColor',Col_100,'LineWidth',2);
s19=scatter(x_hnRNPA3_200,log(y_hnRNPA3_200),72,'MarkerEdgeColor',Col_200,'LineWidth',2);
s20=scatter(x_hnRNPA3_600,log(y_hnRNPA3_600),72,'MarkerEdgeColor',Col_600,'LineWidth',2);

plot(x_vec,log(y_hnRNPA3_0125_fit),'--','Color',Col_0125,'LineWidth',2,'DisplayName',['0.125 {\mu}M']);
plot(x_vec,log(y_hnRNPA3_025_fit),'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,log(y_hnRNPA3_050_fit),'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.5 {\mu}M']);
plot(x_vec,log(y_hnRNPA3_100_fit),'--','Color',Col_100,'LineWidth',2,'DisplayName',['1 {\mu}M']);
plot(x_vec,log(y_hnRNPA3_200_fit),'--','Color',Col_200,'LineWidth',2,'DisplayName',['2 {\mu}M']);
plot(x_vec,log(y_hnRNPA3_600_fit),'--','Color',Col_600,'LineWidth',2,'DisplayName',['6 {\mu}M']);

set(get(get(s15,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s16,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s17,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s18,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s19,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s20,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

 
end

%EWSR1
if EWSR1
s16=scatter(x_EWSR1_025,log(y_EWSR1_025),72,'MarkerEdgeColor',Col_025,'LineWidth',2);
s17=scatter(x_EWSR1_050,log(y_EWSR1_050),72,'MarkerEdgeColor',Col_050,'LineWidth',2);
s18=scatter(x_EWSR1_100,log(y_EWSR1_100),72,'MarkerEdgeColor',Col_100,'LineWidth',2);
s19=scatter(x_EWSR1_200,log(y_EWSR1_200),72,'MarkerEdgeColor',Col_200,'LineWidth',2);

plot(x_vec,log(y_EWSR1_025_fit),'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,log(y_EWSR1_050_fit),'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.5 {\mu}M']);
plot(x_vec,log(y_EWSR1_100_fit),'--','Color',Col_100,'LineWidth',2,'DisplayName',['1 {\mu}M']);
plot(x_vec,log(y_EWSR1_200_fit),'--','Color',Col_200,'LineWidth',2,'DisplayName',['2 {\mu}M']);

set(get(get(s16,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s17,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s18,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s19,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
end

%TAF15
if TAF15
s15=scatter(x_TAF15_0125,log(y_TAF15_0125),72,'MarkerEdgeColor',Col_0125,'LineWidth',2);
s16=scatter(x_TAF15_025,log(y_TAF15_025),72,'MarkerEdgeColor',Col_025,'LineWidth',2);
s17=scatter(x_TAF15_050,log(y_TAF15_050),72,'MarkerEdgeColor',Col_050,'LineWidth',2);
s18=scatter(x_TAF15_100,log(y_TAF15_100),72,'MarkerEdgeColor',Col_100,'LineWidth',2);
s19=scatter(x_TAF15_200,log(y_TAF15_200),72,'MarkerEdgeColor',Col_200,'LineWidth',2);

plot(x_vec,log(y_TAF15_0125_fit),'--','Color',Col_0125,'LineWidth',2,'DisplayName',['0.125 {\mu}M']);
plot(x_vec,log(y_TAF15_025_fit),'--','Color',Col_025,'LineWidth',2,'DisplayName',['0.25 {\mu}M']);
plot(x_vec,log(y_TAF15_050_fit),'--','Color',Col_050,'LineWidth',2,'DisplayName',['0.5 {\mu}M']);
plot(x_vec,log(y_TAF15_100_fit),'--','Color',Col_100,'LineWidth',2,'DisplayName',['1 {\mu}M']);
plot(x_vec,log(y_TAF15_200_fit),'--','Color',Col_200,'LineWidth',2,'DisplayName',['2 {\mu}M']);

set(get(get(s15,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s16,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s17,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s18,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(s19,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

end

p2=plot([0 log(y_min)],[3 log(y_min)],'--k','LineWidth',1); set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');xL=xlabel('R/R^*');
yL=ylabel('log(PR/PR*)');
set(gca,'FontSize',18);
set(gca,'FontWeight','bold');
set(gca,'LineWidth',1);
xL.FontSize=18;
xL.FontWeight='bold';
yL.FontSize=18;
yL.FontWeight='bold';
set(gcf,'color','white');
xlim([x_min x_max]);
ylim([log(y_min) 0])
legend('FontSize',10)
