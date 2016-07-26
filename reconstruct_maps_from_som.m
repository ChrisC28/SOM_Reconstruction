%========================================================================================%
%
% Self Organizing Map Deep Flow Reconstruction
% 
% This script uses a trained map to reconstruct the deep flow using surface observations
% of velocity, ADT and SST (the later strictly not required).
% 
% The script uses gridded products from both AVISO (for surface currents and ADT) and NOAA
% for the SST
% - http://www.aviso.altimetry.fr/fr/donnees/produits/produits-hauteur-de-mer/global.html
% - https://www.ncdc.noaa.gov/oisst
% 
% It also requires a pre-trained map and correlation matrix. These can be obtained by
% running the code som_reconstruction_training.m
% 
% DEPANDANCIES
% The method uses the SOM Toolbox by Esa Alhoniemi,Johan Himberg, Jukka Parviainen and Juha Vesanto.
% http://www.cis.hut.fi/somtoolbox/
% This is free software availble under a GPL licence
% You will need to install this toolbox and add it to the matlab path to use this software
%
% METHOD
% The implementation is identical to that of som_reconstruction_training.m. However, in 
% place of reconstructing the deep flow using the training data, we instead use the 
% maps of U and V over the entire Southern Ocean as inputs. 
% The results of this program are included as Fig. 5 in Chapman and Charontonis (2016)
%
% CONTACT
% This code was written by Chris Chapman and Anatase Charontonis. 
% Questions, commments and bugs can be sent to: 
% chris.chapman.28@gmail.com
% aaclod@locean-ipsl.upmc.fr
% 
%========================================================================================%
%LICENCE 
%Copyright (C) 26/7/2016 by Christopher Chapman and Anastase Alexandre Charontonis 
%
%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
%to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
%OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE 
%USE OR OTHER DEALINGS IN THE SOFTWARE.
%========================================================================================%%  
%

YEAR_TO_RECONSTRUCT = 2007;


%========================================================================================%
%Load the pretrained map file and local correlation matrix
%========================================================================================%
map_file_name = 'som_map_trained_3000.mat'; %file name for later 
load(map_file_name,'sMapobs');

cor_mat_file_name = 'Correlation_Matrix_3000.mat');
load(cor_mat_file_name,'COR_mat');

%========================================================================================%
%LOAD THE SURFACE DATA 
%In this case the data is found the NetCDF files below
%NOTE THE PATHS
%========================================================================================%

aviso_uv_data_path = './';                 
aviso_file_stem = 'uv_dt_SouthernOcean_'; %Surface velocities       

aviso_ssh_data_path = './';
aviso_ssh_file_stem = 'sla_dt_SouthernOcean_'; %Surface SLA

sst_data_path = './';
sst_file_stem = 'sst_SouthernOcean_';           %Sea Surface Temperature

aviso_mean_data_file_name = 'mdt_cnes_cls2013_global.nc';  %Mean Dynamic Topography


%==================================================================================%
% Read in the data to set up the inputs
%==================================================================================%
altimetry_uv_file_to_get = strcat(aviso_file_stem,strcat(num2str(YEAR_TO_RECONSTRUCT),'.nc'));
ssh_file_to_get          = strcat(aviso_ssh_file_stem,strcat(num2str(YEAR_TO_RECONSTRUCT),'.nc'));


lat_mean = ncread(strcat('/net/argos/data/peps/cchlod/AVISO/',aviso_mean_data_file_name),'lat')';  %lat of the MDT may be different to the SLA
lat_grid_sla = ncread(strcat(aviso_uv_data_path,altimetry_uv_file_to_get),'lat')';
lon_grid_sla = ncread(strcat(aviso_uv_data_path,altimetry_uv_file_to_get),'lon')';
time_sla     = ncread(strcat(aviso_uv_data_path,altimetry_uv_file_to_get),'time')';


lat_south_index = find(lat_mean>=lat_grid_sla(1),1);    %Ensure that the MDT and SLA/U and V are on the same grid
lat_north_index = find(lat_mean>=lat_grid_sla(end),1);
lat_count = lat_north_index-lat_south_index+1;


sla_altimetry   = 1.0e-4*ncread(strcat(aviso_ssh_data_path,ssh_file_to_get),'sla');   %Get the SLA. NOTE THE SCALE FACTOR. !!!!
u_sla = ncread(strcat(aviso_uv_data_path,altimetry_uv_file_to_get),'u');              
v_sla = ncread(strcat(aviso_uv_data_path,altimetry_uv_file_to_get),'v');

ssh_mean = squeeze(ncread(strcat('/net/argos/data/peps/cchlod/AVISO/',aviso_mean_data_file_name),'mdt',[1,lat_south_index,1],[Inf,lat_count,Inf],[1,1,1])); %MEAN SSH
u_mean = squeeze(ncread(strcat('/net/argos/data/peps/cchlod/AVISO/',aviso_mean_data_file_name),'u',[1,lat_south_index,1],[Inf,lat_count,Inf],[1,1,1]));     %Mean U
v_mean = squeeze(ncread(strcat('/net/argos/data/peps/cchlod/AVISO/',aviso_mean_data_file_name),'v',[1,lat_south_index,1],[Inf,lat_count,Inf],[1,1,1]));     %Mean V


%Now, get the SST
sst_file_to_get = strcat(sst_file_stem,strcat(int2str(YEAR_TO_RECONSTRUCT),'.nc'));

sst_satellite = 1.0e-2*squeeze(ncread(strcat(sst_data_path,sst_file_to_get),'sst'));
%=============================================================================================%
%Data loaded

addpath('/usr/home/cchlod/somtoolbox')


%Set up the output files
output_file_path = './'
output_file_name = strcat('uv_1000m_from_som_', num2str(YEAR_TO_RECONSTRUCT) ,'.nc');

%Create the file if it doesn't exist
if exist(strcat(output_file_path,output_file_name), 'file') != 2
   nccreate(strcat(output_file_path,output_file_name),'u','Dimensions',{'lon',length(lon_grid_sla),'lat',length(lat_grid_sla),'time',length(time_sla)},'Format','netcdf4_classic')
   nccreate(strcat(output_file_path,output_file_name),'v','Dimensions',{'lon',length(lon_grid_sla),'lat',length(lat_grid_sla),'time',length(time_sla)},'Format','netcdf4_classic')

   nccreate(strcat(output_file_path,output_file_name),'lat','Dimensions',{'lat',length(lat_grid_sla)});
   nccreate(strcat(output_file_path,output_file_name),'lon','Dimensions',{'lon',length(lon_grid_sla)});
   nccreate(strcat(output_file_path,output_file_name),'time','Dimensions',{'time',length(time_sla)});


    ncwrite(strcat(output_file_path,output_file_name),'lat',lat_grid_sla);
    ncwrite(strcat(output_file_path,output_file_name),'lon',lon_grid_sla);
end %if exist(Name, 'file') == 2

for time_step_to_get=1:length(time_sla)

   %=============================================================================================%
   %Set up the data matrix
   %=============================================================================================%

    disp(['Reconstructing flow at time step: ',num2str(time_step_to_get),' of ',num2str(length(time_sla))]);
    [lon_points_sla,lat_points_sla] = meshgrid(lon_grid_sla,lat_grid_sla);

    lon_points_sla = lon_points_sla';
    lat_points_sla = lat_points_sla';

    [n_lat_sla,n_lon_sla] = size(lon_points_sla);


    u_points_sla   = reshape(squeeze(u_sla(:,:,time_step_to_get))+u_mean,[n_lat_sla*n_lon_sla,1]);
    v_points_sla   = reshape(squeeze(v_sla(:,:,time_step_to_get))+v_mean,[n_lat_sla*n_lon_sla,1]);

    adt_points_sla   = reshape(squeeze(sla_altimetry(:,:,time_step_to_get))+ssh_mean,[n_lat_sla*n_lon_sla,1]);
    sst_points_sla   = reshape(squeeze(sst_satellite(:,:,time_step_to_get)),[n_lat_sla*n_lon_sla,1]);

    lon_points_sla = reshape(lon_points_sla,[n_lat_sla*n_lon_sla,1]);
    lat_points_sla = reshape(lat_points_sla,[n_lat_sla*n_lon_sla,1]);

    %Data is packed into a data matrix with the same form as the initial training data. Note the missing deep data in columns 3 and 4
    predictor_data_matrix = [u_points_sla,v_points_sla,zeros([n_lat_sla*n_lon_sla,1]),zeros([n_lat_sla*n_lon_sla,1]),cos(pi/180.0*lon_points_sla ),sin(pi/180.0*lon_points_sla),lat_points_sla,adt_points_sla,sst_points_sla  ];

    %===============================================================================================================%
    %Set up a new map, this time initialised with the test data, we the same topology as the previously trained map

    %extract the hits and the best matching units
    Snew=som_data_struct(predictor_data_matrix(:,:));
    Snew_norm=som_normalize(Snew,sMapobs.comp_norm);

    hits_new=som_hits(sMapobs,Snew_norm);
    bmus_new=som_bmus(sMapobs,Snew_norm);


    %Deal with the empty data
    for i =1:size(Snew_norm.data,1)
        Snew_norm.data(i,isnan(Snew_norm.data(i,:)))=sMapobs.codebook(bmus_new(i),isnan(Snew_norm.data(i,:)));
    end


    %Apply Eqn. 1 of C&C 2016
    for i=1:length(sMapobs.codebook)

        DDist(i,:)= sum(((repmat(sMapobs.codebook(i,:),length(Snew_norm.data),1)-Snew_norm.data).^2).^0.5*diag(COR_mat(i,:)),2)';

    end
    
    %Find the new best-matching units
    for i=1:length(Snew_norm.data)
        wow=find(DDist(:,i)==min(DDist(:,i)));
        bmus_new(i)=wow(1);
    end

    %Reconstruct the data
    sMap2 = sMapobs;
    sMap2 = som_denormalize(sMap2,Sobsn);


    for i =1:size(Snew.data,1)
        Snew.data(i,[3,4])=sMap2.codebook(bmus_new(i),[3,4]);
    end

    %% Grid the data as before 

    u_predicted = double(Snew.data(:,3));
    v_predicted = double(Snew.data(:,4));


    u_predicted = reshape(u_predicted,[n_lat_sla,n_lon_sla]);
    v_predicted = reshape(v_predicted,[n_lat_sla,n_lon_sla]);


    %Write the output
    ncwrite(strcat(output_file_path,output_file_name),'u',u_predicted,[1,1,time_step_to_get]);
    ncwrite(strcat(output_file_path,output_file_name),'v',v_predicted,[1,1,time_step_to_get]);

end %for time_step_to_get

%==============================================================================================%
%END SCRIPT
