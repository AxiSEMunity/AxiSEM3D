%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part of example for AxiSEM3D release paper.
%% This example shows how to create more complex input models for
%% AxiSEM3D than just incorporating uniform anisotropy.
%% Author of this code: Jonathan Wolf, Yale U.
%%
%% In this code we use functions from: 
%% AM Walker and JM Wookey, MSAT - a new toolkit for the analysis of 
%% elastic and seismic anisotropy, 2012.
%%
%% and elastic tensors from:
%% S-i Karato, Deformation of Earth Materials, 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

% multiply strength of US32 by this factor
ani_mult_fac = 3;

%% read elastic tensors and prepare
% elastic tensors from Def of Earth Materials, Karato
C_a = [[236.3 84.5  81.5  0.4   3.4   0.3];
        [84.5   218.5 82.9  -1.8  1.2	0.3];
		[81.5   82.9  208.0 -1.3  6.1   0.2];
		[0.4    -1.8  -1.3  64.9  -0.1  -1.9];
		[3.4    1.2   6.1   -0.1  68.7  0.3];
		[0.3    0.3   0.2   -1.9  0.3   66.6]];
     
C_c =[[223.2   83.6   83.3   0.3   -3.7   0.3];
	     [83.6     209.8  81.9   0.8   1.5    0.3];
		 [83.3     81.9   228.5  0.4   -5.9 0.2];
	     [0.3      0.8    0.4  67.9   0.2   -1.9];
		 [-3.7     1.5	  -5.9	0.2   71.1	0.3];
		 [0.3      0.3    0.2	-1.9  0.3	66.6]];
     
C_e = [[236.8 82.3 84.1 -0.6 0.4 0.1];
	   [82.3 207.7 82.7 -2.6 -0.3 -1.0];
	   [84.1 82.7  217.4  -2.1  -1.9  -0.7];
	   [-0.6  -2.6 -2.1   65.0	0.1	0.4];
       [0.4 -0.3  -1.9 0.1 71.1	-1.4];
       [0.1 -1.0  -0.7 0.4 -1.4 68.5]];
      

%Which Olivine tensor would you like to use?
C = C_e;
rh = 3355;

%plot it for some plausible density
%C = MS_rot3(C,0,0,45);

MS_plot(C, rh, 'quiet', 'polsize', .18, .16, 2.0, 1.0);

%now read US32
M = readtable("US32.csv");
dvs = table2array(M(:,4));
fast_dir = table2array(M(:,5));
magnitude = table2array(M(:,6));

%now read STW105
M = readtable("STW105.csv");
depth_stw05 = table2array(M(:,1))/1000;
density_stw05 = table2array(M(:,2))/1000;
depth_stw05 = 6371.0-flip(depth_stw05);
density_stw05 = flip(density_stw05);

%define model space, see US32.csv file
lon = [-140:1:-50];
lat = [0:1:60];
depth = [10:10:1000];

%interpolate rho to model depths
rho = interp1(depth_stw05,density_stw05,depth);

%% A little crunky but works: Create elastic tensor for every lat-lon-depth
stop_at_depth = 410; %km
stop_index = stop_at_depth/10;

rho_tmp = ones(length(lon), length(lat), stop_index);

C11_tmp = zeros(length(lon), length(lat), stop_index);
C12_tmp = C11_tmp;
C13_tmp = C11_tmp;
C14_tmp = C11_tmp;
C15_tmp = C11_tmp;
C16_tmp = C11_tmp;

C22_tmp = C11_tmp;
C23_tmp = C11_tmp;
C24_tmp = C11_tmp;
C25_tmp = C11_tmp;
C26_tmp = C11_tmp;

C33_tmp = C11_tmp;
C34_tmp = C11_tmp;
C35_tmp = C11_tmp;
C36_tmp = C11_tmp;

C44_tmp = C11_tmp;
C45_tmp = C11_tmp;
C46_tmp = C11_tmp;

C55_tmp = C11_tmp;
C56_tmp = C11_tmp;

C66_tmp = C11_tmp;


depthcount = 0;
%for d = 1:length(depth)
for d = 1:stop_index %only upper mantle
    d
%     ani_perc1 = [];
    latcount = 0;
    for la = 1:length(lat)
    %for la = 1
        
        loncount = 0;
        for lo = 1:length(lon)
             index = (latcount*(length(lon))+lo)+depthcount*(length(lat)*length(lon));
             
             %Calculate fast polarization direction, and ani percentage for vertical
             %raypath. Rotate slightly yo adjust for fast polarization direction
             [pol, perc_ani, ~, ~,~] = MS_phasevels(C,rho(d) * 1000, 90, 0);
             C_loop = MS_rot3(C, 0, 0, -pol);
             
             %calculate isotropic equivalent
             [Ciso] = MS_decomp(C_loop);
             
             %percentage of ani tensor
             int_perc = ani_mult_fac * magnitude(index)/perc_ani;
             
             %make tensor with correct ani strength
             [C_loop, ~] = MS_VRH([int_perc, 1-int_perc], C_loop, rho(d) * 1000, Ciso, rho(d) * 1000);
             
             %adjust fast direction
             C_loop = MS_rot3(C_loop,0,0,fast_dir(index));
             %MS_plot(C_loop, rho(d) * 1000, 'quiet', 'polsize', .18, .16, 2.0, 1.0)
             
             if (MS_checkC(C_loop)==0)
                 fprintf('C corrupted \n')
             end
             
             %populate tensor components
             C11_tmp(lo,la,d) = C_loop(1,1);
             C12_tmp(lo,la,d) = C_loop(1,2);
             C13_tmp(lo,la,d) = C_loop(1,3);
             C14_tmp(lo,la,d) = C_loop(1,4);
             C15_tmp(lo,la,d) = C_loop(1,5);
             C16_tmp(lo,la,d) = C_loop(1,6);

             C22_tmp(lo,la,d) = C_loop(2,2);
             C23_tmp(lo,la,d) = C_loop(2,3);
             C24_tmp(lo,la,d) = C_loop(2,4);
             C25_tmp(lo,la,d) = C_loop(2,5);
             C26_tmp(lo,la,d) = C_loop(2,6);             
             
             C33_tmp(lo,la,d) = C_loop(3,3);
             C34_tmp(lo,la,d) = C_loop(3,4);
             C35_tmp(lo,la,d) = C_loop(3,5);
             C36_tmp(lo,la,d) = C_loop(3,6);      
             
             C44_tmp(lo,la,d) = C_loop(4,4);
             C45_tmp(lo,la,d) = C_loop(4,5);
             C46_tmp(lo,la,d) = C_loop(4,6);        
             
             C55_tmp(lo,la,d) = C_loop(5,5);
             C56_tmp(lo,la,d) = C_loop(5,6);                
             
             C66_tmp(lo,la,d) = C_loop(6,6);                
             
             rho_tmp(lo,la,d) = rho(d) * 1000;
             loncount = loncount+1;
        end
        latcount = latcount+1;
    end
    depthcount = depthcount+1;
%     ani_perc1_depth = [ani_perc1_depth, mean(ani_perc1)];
end



%% Write NETCDF4
ncid = netcdf.create('azi_410_3.nc', 'NETCDF4');
dimdep = netcdf.defDim(ncid,'dimdep',stop_index);
dimlat = netcdf.defDim(ncid,'dimlat',length(lat));
dimlon = netcdf.defDim(ncid,'dimlon',length(lon));

lats = netcdf.defVar(ncid,'latitude','NC_FLOAT',dimlat);
netcdf.putVar(ncid,lats,lat);
longs = netcdf.defVar(ncid,'longitude','NC_FLOAT',dimlon);
netcdf.putVar(ncid,longs,lon);
depths = netcdf.defVar(ncid,'depth','NC_FLOAT',dimdep);
netcdf.putVar(ncid,depths,depth(1:stop_index));

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

netcdf.putVar(ncid,c11, C11_tmp);
netcdf.putVar(ncid,c12, C12_tmp);
netcdf.putVar(ncid,c13, C13_tmp);
netcdf.putVar(ncid,c14, C14_tmp);
netcdf.putVar(ncid,c15, C15_tmp);
netcdf.putVar(ncid,c16, C16_tmp);

netcdf.putVar(ncid,c22, C22_tmp);
netcdf.putVar(ncid,c23, C23_tmp);
netcdf.putVar(ncid,c24, C24_tmp);
netcdf.putVar(ncid,c25, C25_tmp);
netcdf.putVar(ncid,c26, C26_tmp);

netcdf.putVar(ncid,c33, C33_tmp);
netcdf.putVar(ncid,c34, C34_tmp);
netcdf.putVar(ncid,c35, C35_tmp);
netcdf.putVar(ncid,c36, C36_tmp);

netcdf.putVar(ncid,c44, C44_tmp);
netcdf.putVar(ncid,c45, C45_tmp);
netcdf.putVar(ncid,c46, C46_tmp);

netcdf.putVar(ncid,c55, C55_tmp);
netcdf.putVar(ncid,c56, C56_tmp);

netcdf.putVar(ncid,c66, C66_tmp);
% close the file
netcdf.close(ncid);