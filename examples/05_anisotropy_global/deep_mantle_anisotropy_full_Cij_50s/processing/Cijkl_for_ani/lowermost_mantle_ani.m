%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code is an example for how to create AxiSEM3D netcdf input for an 
%% elastic tensor Cij. 
%% In the following, we incorporate the anisotropy from anisotropic PREM 
%% in the upper mantle.
%%
%% For questions, please contact Jonathan Wolf (jonathan.wolf@yale.edu)
%%
%% Jonathan Wolf and Andrea Tesoniero have contributed to this code.
%% There is no guarantee for functionality or fitness for any particular
%% purpose.
%%
%% In this code we use functions from: 
%% AM Walker and JM Wookey, MSAT - a new toolkit for the analysis of 
%% elastic and seismic anisotropy, 2012.
%%
%% and elastic tensors from:
%% NM Creasy and L Miyagi and MD Long, A Library of Elastic Tensors for 
%% Lowermost Mantle Seismic Anisotropy Studies and Comparison With Seismic 
%% Observations, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%specify depth and lon, lat
depth = [2641, 2890];
lon = 0:30:180;
lat = -90:30:90;

%The elastic tensors are all from the Creasy et al. (2020) publication

%Wookey et al., 2005b 126 Gpa, 2800 K				Density (kg/m^3): 5191
%Table S1
C_br = [[1013.272	448.399	450.036	-0.071	-0.129	18.892];
[448.399	933.315	445.996	-0.675	-0.832	23.349];
[450.036	445.996	982.483	-0.471	-0.439	1.670];
[-0.071	-0.675	-0.471	256.539	9.455	0.258];
[-0.129	-0.832	-0.439	9.455	274.783	-0.102];
[18.892	23.349	1.670	0.258	-0.102	262.135]];
C_br = MS_rot3(C_br, 0, 0, 0);
rh_br = 5191; %Density of bridgmanite

%Karki et al., 1999 135 Gpa, 3000 K				Density (kg/m^3): 5035	
%from Table S3
C_fe = [[1003.571	414.375	419.063	-0.893	0.418	40.209];
[414.375	1011.121	411.513	1.219	-0.177	-46.414];
[419.063	411.513	1006.433	-0.325	-0.242	6.206];
[-0.893	1.219	-0.325	272.222	5.706	-0.340];
[0.418	-0.177	-0.242	5.706	278.409	-0.844];
[40.209	-46.414	6.206	-0.340	-0.844	280.716]];
C_fe = MS_rot3(C_fe, 0, 0, 0);
rh_fe = 5035;

%Stackhouse et al., 2005 136 GPa, 3000 K	Density (kg/m^3): 5336	
%from Table S7
C_ppv = [[1041.977	447.453	425.724	0.273	-0.395	5.270];
[447.453	1041.731	449.124	0.784	-0.305	-6.771];
[425.724	449.124	1041.969	-0.454	0.875	-14.462];
[0.273	0.784	-0.454	291.640	9.254	0.459];
[-0.395	-0.305	0.875	9.254	306.521	0.048];
[5.270	-6.771	-14.462	0.459	0.048	302.402]];

C_ppv = MS_rot3(C_ppv, 0, 0, -25);
C_ppv = MS_rot3(C_ppv, 0, -50, 0);
C_ppv = MS_rot3(C_ppv, 0, 80, 0);
rh_ppv = 5336;


%Wookey et al., 2005b 126 Gpa, 2800 K; Density (kg/m^3): 5152	
%TableS9
C_brfp = [[1006.982	440.176	445.370	-0.485	0.199	8.503];
[440.176	949.262	439.990	0.060	-0.494	21.353];
[445.370	439.990	982.924	-0.517	-0.781	2.990];
[-0.485	0.060	-0.517	262.339	9.091	0.342];
[0.199	-0.494	-0.781	9.091	278.861	-0.510];
[8.503	21.353	2.990	0.342	-0.510	265.341]];
C_brfp = MS_rot3(C_brfp, 0, 0, 0);
rh_brfp = 5152; %Density of Br

%Stackhouse et al., 2005 136 GPa, 3000 K; Density (kg/m^3): 5261	
%from Table S16
C_ppvfp = [[1031.804	439.525	423.810	-0.020	-0.170	14.285];
[439.525	1033.351	439.831	0.885	-0.285	-16.946];
[423.810	439.831	1032.780	-0.414	0.588	-9.263];
[-0.020	0.885	-0.414	286.656	8.317	0.257];
[-0.170	-0.285	0.588	8.317	299.234	-0.192];
[14.285	-16.946	-9.263	0.257	-0.192	296.276]];
C_ppvfp = MS_rot3(C_ppvfp, 0, 0, 0); %this one won
rh_ppvfp = 5261;

%Plot the tensors
MS_plot(C_br, rh_br, 'quiet', 'polsize', .18, .16, 2.0, 1.0)
MS_plot(C_fe, rh_fe, 'quiet', 'polsize', .18, .16, 2.0, 1.0)
MS_plot(C_ppv, rh_ppv, 'quiet', 'polsize', .18, .16, 2.0, 1.0)
MS_plot(C_br, rh_brfp, 'quiet', 'polsize', .18, .16, 2.0, 1.0)
MS_plot(C_fe, rh_ppvfp, 'quiet', 'polsize', .18, .16, 2.0, 1.0)

C = NaN(6,6,length(lat) * length(lon), length(depth));

%choose elastic tensor here. We use a Ppv tensor here.
for i=1:length(depth)
    for j = 1:length(lat) * length(lon)        
        C(:,:,j,i) = C_ppv(:,:);        
    end
end

Cijkl = C;

%% Write NETCDF4
ncid = netcdf.create('lowermost_mantle_ani.nc', 'NETCDF4');
dimdep = netcdf.defDim(ncid,'dimdep',length(depth));
dimlat = netcdf.defDim(ncid,'dimlat',length(lat));
dimlon = netcdf.defDim(ncid,'dimlon',length(lon));

lats = netcdf.defVar(ncid,'latitude','NC_FLOAT',dimlat);
netcdf.putVar(ncid,lats,lat);
longs = netcdf.defVar(ncid,'longitude','NC_FLOAT',dimlon);
netcdf.putVar(ncid,longs,lon);
depths = netcdf.defVar(ncid,'depth','NC_FLOAT',dimdep);
netcdf.putVar(ncid,depths,depth);

c11 = netcdf.defVar(ncid, 'C11','NC_FLOAT',[dimlon,dimlat,dimdep]);
c12 = netcdf.defVar(ncid, 'C12','NC_FLOAT',[dimlon,dimlat,dimdep]);
c13 = netcdf.defVar(ncid, 'C13','NC_FLOAT',[dimlon,dimlat,dimdep]);
c14 = netcdf.defVar(ncid, 'C14','NC_FLOAT',[dimlon,dimlat,dimdep]);
c15 = netcdf.defVar(ncid, 'C15','NC_FLOAT',[dimlon,dimlat,dimdep]);
c16 = netcdf.defVar(ncid, 'C16','NC_FLOAT',[dimlon,dimlat,dimdep]);

c22 = netcdf.defVar(ncid, 'C22','NC_FLOAT',[dimlon,dimlat,dimdep]);
c23 = netcdf.defVar(ncid, 'C23','NC_FLOAT',[dimlon,dimlat,dimdep]);
c24 = netcdf.defVar(ncid, 'C24','NC_FLOAT',[dimlon,dimlat,dimdep]);
c25 = netcdf.defVar(ncid, 'C25','NC_FLOAT',[dimlon,dimlat,dimdep]);
c26 = netcdf.defVar(ncid, 'C26','NC_FLOAT',[dimlon,dimlat,dimdep]);

c33 = netcdf.defVar(ncid, 'C33','NC_FLOAT',[dimlon,dimlat,dimdep]);
c34 = netcdf.defVar(ncid, 'C34','NC_FLOAT',[dimlon,dimlat,dimdep]);
c35 = netcdf.defVar(ncid, 'C35','NC_FLOAT',[dimlon,dimlat,dimdep]);
c36 = netcdf.defVar(ncid, 'C36','NC_FLOAT',[dimlon,dimlat,dimdep]);

c44 = netcdf.defVar(ncid, 'C44','NC_FLOAT',[dimlon,dimlat,dimdep]);
c45 = netcdf.defVar(ncid, 'C45','NC_FLOAT',[dimlon,dimlat,dimdep]);
c46 = netcdf.defVar(ncid, 'C46','NC_FLOAT',[dimlon,dimlat,dimdep]);

c55 = netcdf.defVar(ncid, 'C55','NC_FLOAT',[dimlon,dimlat,dimdep]);
c56 = netcdf.defVar(ncid, 'C56','NC_FLOAT',[dimlon,dimlat,dimdep]);

c66 = netcdf.defVar(ncid, 'C66','NC_FLOAT',[dimlon,dimlat,dimdep]);

c11_tmp = ones(length(lon), length(lat), length(depth));
c12_tmp = c11_tmp;
c13_tmp = c11_tmp;
c14_tmp = c11_tmp;
c15_tmp = c11_tmp;
c16_tmp = c11_tmp;

c22_tmp = c11_tmp;
c23_tmp = c11_tmp;
c24_tmp = c11_tmp;
c25_tmp = c11_tmp;
c26_tmp = c11_tmp;

c33_tmp = c11_tmp;
c34_tmp = c11_tmp;
c35_tmp = c11_tmp;
c36_tmp = c11_tmp;

c44_tmp = c11_tmp;
c45_tmp = c11_tmp;
c46_tmp = c11_tmp;
c55_tmp = c11_tmp;
c56_tmp = c11_tmp;

c66_tmp = c11_tmp;

%populate the variables (a little cruncky but works just fine)

for i = 1:length(depth)
    
    c11_tmp(:,:,i) = reshape(Cijkl(1,1,:,i),[length(lat), length(lon)])' ;
    
    c12_tmp(:,:,i) = reshape(Cijkl(1,2,:,i),[length(lat), length(lon)])' ;
    c13_tmp(:,:,i) = reshape(Cijkl(1,3,:,i),[length(lat), length(lon)])' ;
    c14_tmp(:,:,i) = reshape(Cijkl(1,4,:,i),[length(lat), length(lon)])' ;
    c15_tmp(:,:,i) = reshape(Cijkl(1,5,:,i),[length(lat), length(lon)])' ;
    c16_tmp(:,:,i) = reshape(Cijkl(1,6,:,i),[length(lat), length(lon)])' ;
    
    c22_tmp(:,:,i) = reshape(Cijkl(2,2,:,i),[length(lat), length(lon)])' ;
    c23_tmp(:,:,i) = reshape(Cijkl(2,3,:,i),[length(lat), length(lon)])' ;
    c24_tmp(:,:,i) = reshape(Cijkl(2,4,:,i),[length(lat), length(lon)])' ;
    c25_tmp(:,:,i) = reshape(Cijkl(2,5,:,i),[length(lat), length(lon)])' ;
    c26_tmp(:,:,i) = reshape(Cijkl(2,6,:,i),[length(lat), length(lon)])' ;
    
    c33_tmp(:,:,i) = reshape(Cijkl(3,3,:,i),[length(lat), length(lon)])' ;
    c34_tmp(:,:,i) = reshape(Cijkl(3,4,:,i),[length(lat), length(lon)])' ;
    c35_tmp(:,:,i) = reshape(Cijkl(3,5,:,i),[length(lat), length(lon)])' ;
    c36_tmp(:,:,i) = reshape(Cijkl(3,6,:,i),[length(lat), length(lon)])' ;
    
    c44_tmp(:,:,i) = reshape(Cijkl(4,4,:,i),[length(lat), length(lon)])' ;
    c45_tmp(:,:,i) = reshape(Cijkl(4,5,:,i),[length(lat), length(lon)])' ;
    c46_tmp(:,:,i) = reshape(Cijkl(4,6,:,i),[length(lat), length(lon)])' ;
    
    c55_tmp(:,:,i) = reshape(Cijkl(5,5,:,i),[length(lat), length(lon)])' ;
    c56_tmp(:,:,i) = reshape(Cijkl(5,6,:,i),[length(lat), length(lon)])' ;
    
    c66_tmp(:,:,i) = reshape(Cijkl(6,6,:,i),[length(lat), length(lon)])' ;
    
end


netcdf.putVar(ncid,c11, c11_tmp);

netcdf.putVar(ncid,c12, c12_tmp);
netcdf.putVar(ncid,c13, c13_tmp);
netcdf.putVar(ncid,c14, c14_tmp);
netcdf.putVar(ncid,c15, c15_tmp);
netcdf.putVar(ncid,c16, c16_tmp);

netcdf.putVar(ncid,c22, c22_tmp);
netcdf.putVar(ncid,c23, c23_tmp);
netcdf.putVar(ncid,c24, c24_tmp);
netcdf.putVar(ncid,c25, c25_tmp);
netcdf.putVar(ncid,c26, c26_tmp);

netcdf.putVar(ncid,c33, c33_tmp);
netcdf.putVar(ncid,c34, c34_tmp);
netcdf.putVar(ncid,c35, c35_tmp);
netcdf.putVar(ncid,c36, c36_tmp);

netcdf.putVar(ncid,c44, c44_tmp);
netcdf.putVar(ncid,c45, c45_tmp);
netcdf.putVar(ncid,c46, c46_tmp);

netcdf.putVar(ncid,c55, c55_tmp);
netcdf.putVar(ncid,c56, c56_tmp);

netcdf.putVar(ncid,c66, c66_tmp);
% close the file
netcdf.close(ncid);
%ncdisp('prem_ani_MATLAB.nc')
