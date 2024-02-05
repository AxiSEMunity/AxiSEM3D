%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code is an example for how to create AxiSEM3D netcdf input for an 
%% elastic tensor Cij. 
%% In the following, we incorporate the anisotropy from anisotropic PREM 
%% in the upper mantle.
%%
%% For questions, please contact Jonathan Wolf (jonathan.wolf@yale.edu)
%%
%% Andrea Tesoniero and Jonathan Wolf have contributed to this code.
%% There is no guarantee for functionality or fitness for any particular
%% purpose.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% read in PREM
prem = importdata('PREM_1s.csv' ,',');
% radius,depth,density,Vpv,Vph,Vsv,Vsh,eta,Q-mu,Q-kappa
Cijkl = zeros(6, 6, length(prem));
depth = NaN(1,length(prem));
rh = NaN(1,length(prem));
%only read the PREM values for the layer that you are interested in
n = 1;  
for i = 1:length(prem)
    if (prem(i,8) ~= 1.0) % eta
        if prem(i,6) > 0. %  not the outer core
           
           % These are the so-called Love parameters.
           A = prem(i,3) *  prem(i,5) * prem(i,5); % A = rho * vpv^2
           C = prem(i,3) *  prem(i,4) * prem(i,4); % C = rho * vph^2
           L = prem(i,3) *  prem(i,6) * prem(i,6); % L = rho * vsv^2
           N = prem(i,3) *  prem(i,7) * prem(i,7); % N = rho * vsh^2
           F = prem(i,8) * (A - 2 * L);
           
           % Translate Love parameters to Cij
           Cijkl(1,1,n) = A; Cijkl(2,2,n) = A;
           Cijkl(1,2,n) =  A - 2 * N; Cijkl(2,1,n) = A - 2 * N;
           Cijkl(1,3,n) = F; Cijkl(2,3,n) = F; Cijkl(3,1,n) = F; Cijkl(3,2,n) = F;
           Cijkl(3,3,n) = C;
           Cijkl(4,4,n) = L; Cijkl(5,5,n) = L;
           Cijkl(6,6,n) = N;
           depth(n) = prem(i,2);
           rh(n) = prem(i,3);
           vs(n) = prem(i,6); %just in case this is of interest. We don't use it here
           n = n + 1;
        end
    end
end

% trim Cijkls and depth
for i = 1:length(prem)
    if ~any(Cijkl(:,:,i))
        break
    end
end

% make sure, Cijkl is only defined for the depth range in which PREM is
% anisotropic
Cijkl(:,:,i:end)= [];
depth(:,i:end) = [];

% AxiSEM3D does not like it if you specify the same depth twice, which is
% what PREM does at discontinuities. Avoid problems by separating layers.
for i=2:length(depth)
    if depth(i-1) == depth(i)
        depth(i)=depth(i)+0.1;
    end
end

% Specify whether you would like to include anisotropic PREM for only
% regions of Earth.
% Edit lon and lat as you wish!
% Of course, this looks strange if we are doing the whole globe anyways
lon = -180:30:180;
lat =  -90:30:90;

C = NaN(6,6,length(lat) * length(lon), length(depth));


for i=1:length(depth)
    for j = 1:length(lat) * length(lon)        
        C(:,:,j,i) = Cijkl(:,:,i);        
    end
end

Cijkl = C;


%% Write NETCDF4
ncid = netcdf.create('prem_ani_upper.nc', 'NETCDF4');
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

%This ugly, but it works.

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
