%========================================================================================%
%
% Self Organizing Map Reconstruction
%
% This script demonstrates the implementation of the method for the reconstruction of 
% deep velocity from surface observations described in Chapman & Charontonis (2016).
% In particular, this script demonstrates the:
% - training of the SOM;
% - the modification of the SOM based on local correlations (Eqn. 1 in Chapman & Charontonis 2016)
% - the reconstruction of deep velocities and the validation (Fig. 3 in C&C2016)
%
% DEPANDANCIES
% The method uses the SOM Toolbox by Esa Alhoniemi,Johan Himberg, Jukka Parviainen and Juha Vesanto.
% http://www.cis.hut.fi/somtoolbox/
% This is free software availble under a GPL licence
% You will need to install this toolbox and add it to the matlab path to use this software
% 
% Additionally, you will need access to the ANDRO float data set (Ollitraud & Rannou 2013)
% http://www.umr-lops.fr/Donnees/ANDRO
% and atleast some measure of the surface velocity and ADT interpolated to the ANDRO data
% locations and times.  
% In this study we have used AVISO for velocity and ADT: 
% - http://www.aviso.altimetry.fr/fr/donnees/produits/produits-hauteur-de-mer/global.html
% and the NOAA Optimum Interpolation Sea Surface Temperature
% - https://www.ncdc.noaa.gov/oisst
%
% For the purposes of this script, the data interpolated to the 
% deep velocity locations, plus the deep velocity themselves, are stored in the files: 
% adt_interp_argo_sst.mat 
% and 
% oisst_interp_argo_sst.mat
%
% The later is strictly not necessary, but it does improve the reconstruction.
% 
% METHOD
% As described in Chapman & Charontonis, the SOM is trained using current velocity  at
% depth, and ADT, surface currents and SST at the surface from satellites. 
% 20% of the data set is retained for validaion. 
% - The map is trained iteratively, using decreasing left and right influence radii. 
% - We then apply the local correlation weighting IN DATA SPACE!!! (line )   
% - The deep velocity is reconstructed by mapping the surface data to the best-matching units
% of the modified SOM
%
% CONTACT
% This code was written by Chris Chapman and Anatase Charontonis. 
% Questions, commments and bugs can be sent to: 
% chris.chapman.28@gmail.com
% aaclod@locean-ipsl.upmc.fr
% 
%%========================================================================================%
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
%%========================================================================================%

%Do you want figures to appear? 
SHOW_FIGURES = 0

%==============================================================================%
% LOAD THE IMPUT DATA
%load the input data that includes the colocated ssh derived velocities and the 
%deep velocities from the Argo data sets
%==============================================================================%
base_path = '/net/argos/data/peps/cchlod/AVISO/'
file_name = 'adt_interp_argo_ssh.mat'
sst_file_path  = '/net/argos/data/peps/cchlod/OISST/'
sst_file_name  = 'oisst_interp_argo_sst.mat'

%================================================================%
% Load the ADT
input_file_object = load(strcat(base_path,file_name),'-mat','');
pressure          = input_file_object.andro_pressure_parking;
lat_points        = input_file_object.andro_lat;
lon_points        = input_file_object.andro_lon;

deep_u_vel        = input_file_object.u_vel;
deep_v_vel        = input_file_object.v_vel;


deep_u_err        = input_file_object.u_vel_err;
deep_v_err        = input_file_object.v_vel_err;

surf_u_vel        = input_file_object.u_surf_mean_andro + input_file_object.u_surf_anom_andro;
surf_v_vel        = input_file_object.v_surf_mean_andro + input_file_object.v_surf_anom_andro;

sla        = input_file_object.sla_andro;
mdt        = input_file_object.mdt_andro;
adt        = mdt + sla;
dates      = input_file_object.andro_dates_deep;
pressure   = input_file_object.andro_pressure_parking;

%================================================================%
clear input_file_object

% Load the SST
input_file_object = load(strcat(sst_file_path,sst_file_name),'-mat','');
sst               = input_file_object.sst_andro;
clear input_file_object
%================================================================%

%================================================================%
% Fixed parameters for the training of the map
% We've tested numerous values here and the work well, but they
% can be modified without fear of breaking the script
radius_left = [6,4,2,1,0.5,0.1];
radius_right = [4,2,1,0.5,0.1,0.0001];
training_length = [60,100,200,200,140,60];

n_munits = 3000; %number of classes                                
%================================================================%
% if the somtoolbox is not in the path, add it here
addpath('/usr/home/cchlod/somtoolbox')
%form data matrix
rng('default');
rng(10);

%We select 80% of the data set as training data, the remaider is 
%retained for validations
percent_training = 0.8;

%Domain bounds
lon_west = 0.0;
lon_east = 360.0;

lat_south = -70.0; %Typical Southern Ocean lattiudes
lat_north = -35.0;

dates_ymd = datevec(dates);
dates_ymd = dates_ymd(:,1:3);

argo_years        = dates_ymd(:,1); %We convert the dates to YMD format from "days since" formate


%==========================================================================================================================================================================================
% We restrict our attention to observations in the Southern Ocean (for now)
% all points that fall outside of the box defined by lat_north/lat_south/lon_west/lon_east and off a certain pressure surface are rejected
%==========================================================================================================================================================================================

indicies_in_box = find((lat_points>lat_south & lat_points<lat_north) & (lon_points>lon_west & lon_points<lon_east) & (argo_years>2004 & argo_years<2011) & (pressure>950 & pressure<1050 ) );

local_lon = lon_points(indicies_in_box);
local_lat = lat_points(indicies_in_box);
dates_local = dates(indicies_in_box);
local_pressure = pressure(indicies_in_box ); 


dates_ymd = datevec(dates_local);
dates_ymd = dates_ymd(:,1:3);


deep_u_vel_local = deep_u_vel(indicies_in_box);
deep_v_vel_local = deep_v_vel(indicies_in_box);


surf_u_vel_local = surf_u_vel(indicies_in_box);
surf_v_vel_local = surf_v_vel(indicies_in_box);

deep_u_err_local  = deep_u_err(indicies_in_box); 
deep_v_err_local  = deep_v_err(indicies_in_box); 

sla_local = sla(indicies_in_box);
mdt_local = mdt(indicies_in_box);
sst_local = sst(indicies_in_box);

adt_local = sla_local+mdt_local;


years        = dates_ymd;
years(:,2:3) = 0.0;

%===========================================================================================================%
%The map cannot understand spherical geometry and does not think that points at 359 degrees of longitude and 
% 1 degree of longitude are close in space
% Thus we use a poor-man's quaterion transform and give the map the cos and sin of the lonitude  
%===========================================================================================================%
cos_lon = cos((pi/180.0)*local_lon);
sin_lon = sin((pi/180.0)*local_lon);


%Split the data into training and test datasets
n_points_in_region = length(local_lon);
size_of_sample     = floor( n_points_in_region * percent_training);

training_sample_indicies = randsample(n_points_in_region,size_of_sample);
test_sample_indicies     = setdiff([1:n_points_in_region],training_sample_indicies);

%form the training data matrix (N_trianing x 9) 
training_data_matrix = [surf_u_vel_local(training_sample_indicies),surf_v_vel_local(training_sample_indicies),deep_u_vel_local(training_sample_indicies),deep_v_vel_local(training_sample_indicies),...
                        cos_lon(training_sample_indicies),sin_lon(training_sample_indicies),local_lat(training_sample_indicies),adt_local(training_sample_indicies),sst_local(training_sample_indicies)];

%form the test data matrix (N_test x 9) 
test_data_matrix     = [surf_u_vel_local( test_sample_indicies ),surf_v_vel_local(test_sample_indicies),deep_u_vel_local(test_sample_indicies),deep_v_vel_local(test_sample_indicies),...
                        cos_lon(test_sample_indicies),sin_lon(test_sample_indicies),local_lat(test_sample_indicies),adt_local(test_sample_indicies),sst_local(test_sample_indicies)];


[n_traning_points,n_vars] = size(training_data_matrix); 

%=============================================================================================================%
%Here we invoke the SOM
%=============================================================================================================% 

Sobs=som_data_struct(training_data_matrix(:,:)); %Set up the SOM data structure
Sobsn=som_normalize(Sobs,'var');                 %Normalize the data 


map_file_name = strcat('som_map_trained_',num2str(n_munits),'_classes.mat'); %file name for later 

disp(['Training map with ', num2str(n_munits), ' classes' ])

%Make a map object
sMapobs=som_make(Sobsn,'munits',n_munits);

for i = 1:length(radius_left)
    
    %Iteratively train the SOM with progressively longer training lengths and smaller influence radii
    sMapobs=som_batchtrain(sMapobs,Sobsn,'radius',[radius_left(i) radius_right(i)],'trainlen',training_length(i),'tracking',0);

end %for i = 1:length(radius_left) 
   
%Saving the map after training is a good idea
save(map_file_name,'sMapobs');
 

if SHOW_FIGURES 
    fig1 =figure(1)
    som_show(sMapobs) 
end %if SHOW FIGURES

%Get the hits map and the best matching units
hits=som_hits(sMapobs,Sobsn);
bmus=som_bmus(sMapobs,Sobsn);


if SHOW_FIGURES 
    fig2 = figure(2)
    som_cplane(sMapobs,hits)
    colorbar()

end %if SHOW FIGURES

%Set the data matrix to NaN places where the best-matching units are empty
for i =1:size(Sobsn.data,1)
        Sobsn.data(i,isnan(Sobsn.data(i,:)))=sMapobs.codebook(bmus(i),isnan(Sobsn.data(i,:)));
end

%=========================================================================================================%
%LOCAL CORRELATIONS IN DATA SPACE
% Here, we apply Eqn. 1 in C&C 2016 to the trained map. 
%=========================================================================================================% 
COR_mat = 99*ones([length(sMapobs.codebook),n_vars]);
for i =1:length(sMapobs.codebook)
 
    if ~isempty(Sobsn.data(bmus==i,:))
        %Calculate the local correlations in the data space
        corr_train=corr(Sobsn.data(bmus==i,:));
        COR_mat(i,:)=1+sqrt(corr_train(:,3).^2+corr_train(:,4).^2);

        end %if isempty
end %for i=1:length(sMapobs.codebook) 

%Save the correlation matrix
cor_mat_file_name = strcat('Correlation_Matrix_',num2str(n_munits),'.mat');
save(cor_mat_file_name,'COR_mat');
%Set up a new map, this time initialised with the test data, we the same topolog as the previously trained map
Snew=som_data_struct(test_data_matrix(:,:));
Snew_norm=som_normalize(Snew,Sobsn);   

for i =1:size(Snew_norm.data,1)
    Snew_norm.data(i,isnan(Snew_norm.data(i,:)))=sMapobs.codebook(bmus(i),isnan(Snew_norm.data(i,:)));
end %for i =1:size(Snew_norm.data,1)


for i=1:length(sMapobs.codebook)
        %Find the weighted distance in data space
        DDist(i,:)= sum(((repmat(sMapobs.codebook(i,:),length(Snew_norm.data),1)-Snew_norm.data).^2).^0.5*diag(COR_mat(i,:)),2)';    
    
end %i=1:length(sMapobs.codebook)

%The new best matching unit is with the smallest weighted distance in data space for the inputs
for i=1:length(Snew_norm.data)
    wow=find(DDist(:,i)==min(DDist(:,i)));
    bmus_new(i)=wow(1);
end


%==========================================================%
% RECONSTRUCTION
% Finally, we reconstruct the deep observations from the 
% surface training data
%==========================================================%

%Again another data structure
sMap_reconstruct = sMapobs;
%Denormalize back to physical space
sMap_reconstruct = som_denormalize(sMap_reconstruct,Sobsn);

%Use the new best matching units to fill in the "missing" deep velocity
%(colums 3 and 4 in the training and test matricies
%
for i =1:size(Snew.data,1)
    Snew.data(i,[3,4])=sMap_reconstruct.codebook(bmus_new(i),[3,4]);
end
%%

%Voila! 
u_predicted = Snew.data(:,3);
v_predicted = Snew.data(:,4);


%Compute the speed error at each point in training dataset
error_speed = sqrt((u_predicted-deep_u_vel_local(test_sample_indicies)).^2 + (v_predicted-deep_v_vel_local(test_sample_indicies)).^2);

%Compute the total RMSE of the speed error
total_rmse = sqrt( nanmean((u_predicted-deep_u_vel_local(test_sample_indicies)).^2 + (v_predicted-deep_v_vel_local(test_sample_indicies) ).^2) )
disp(['Total Speed RMSE: ', num2str(total_rmse) ] )
if SHOW_FIGURES 
    fig4=figure(4)
    scatter(u_predicted,deep_u_vel_local(test_sample_indicies))
    hold on
    plot(deep_u_vel_local(test_sample_indicies),deep_u_vel_local(test_sample_indicies),'r')
    title('U predicted vs. U observed')
    hold off

    fig5=figure(5)
    scatter(v_predicted,deep_v_vel_local(test_sample_indicies))
    hold on
    plot(deep_v_vel_local(test_sample_indicies),deep_v_vel_local(test_sample_indicies),'r')
    title('V predicted vs. V observed')
    hold off
end %if SHOW_FIGURES
%================================================================================================%
%END OF SCRIPT
%================================================================================================%
