function prior=getMapPrior(lat,lon,G)


%% Main - compute image label prior
% 
% Input: location, heading, camera intrinsics
% Output: projected images (with label area)        
% 
% figures - all map, projected image, image label prior
%
% Jifu Zhou (jifuz@andrew.cmu.edu) 11/11/2012

%close all



%% location / heading
%location
%loc=[40.444524,-79.946276];%CICEntrance
%loc=[40.444604,-79.944001];%morewood
%loc=[40.444628,-79.943631];%morewood
%loc=[40.444636,-79.943186];%morewood
%loc=[40.442513,-79.937977];%a building on MM and Forbes
loc=[lat,lon];
%Gaze direction
%G=[10 0 0];


%% camera parameters
lat=loc(1);
lon=loc(2);
Cam=cameraParam(lat,lon,G);


%% load data 
[osm,outputfilename]=retrieveData(lat,lon,G);


%% specify information structure - visiblility level of each obj
global lowestLevel
global VisLevel
lowestLevel=10;
VisLevel=1;
[osm_s]=osmInfoStruct(osm);


%% visualize complete data
figure
axis equal
hold on
plot(0,0,'x','MarkerSize',50,'LineWidth',3)
visualizeOSM(Cam,osm_s);


%% filter - retrieve ojects in the field of view, in order !
%% compute projection
reprojection=computeProjection(Cam,osm_s);
visualizeProjection(reprojection,osm_s,Cam)


%% compute label area of projected images !
[labelLayers]=visualizeProjection_fillingLabels(reprojection,osm_s,Cam);
[prior]=mergeLayers(labelLayers);
figure
imagesc(prior)

outputfilename=strrep(outputfilename,'.xml','.mat');
outputfilename=strrep(outputfilename,'rawOSM','priors');
save(outputfilename,'prior')



end