% Main - compute image label prior
% 
% Input: location, heading, camera intrinsics
% Output: projected images (with label area)        
% 
% figures - all map, projected image, image label prior
%
% Jifu Zhou (jifuz@andrew.cmu.edu) 11/11/2012

 close all



%% location / heading
%location
%loc=[40.444524,-79.946276];%CICEntrance
%loc=[40.444604,-79.944001];%morewood
%loc=[40.444628,-79.943631];%morewood
%loc=[40.444636,-79.943186];%morewood
%loc=[40.442513,-79.937977];%a building on MM and Forbes
%loc=[40.444465,-79.942025];%Forbes turnaround
% loc=[40.445179,-79.948778];%winthrop /bustop on the right
% loc=[40.445537,-79.948825];%winthrop /CLOSER bustop on the right
%loc=[40.445072,-79.948762];%intersection winthrop
loc=[48.99921658619, 8.4746294522095];%kitti 0
loc=[48.999193021883 8.4744379792903];%kitti 20
loc=[48.999204386904 8.474534643113];%kitti 10
%loc=[40.44461,-79.943943];%morewood before parking lot

t=(-2.9185069803847/pi)*180;
G1=rotateVectorClockwise([1 0],-t);
G=[G1' 0];
%Gaze direction
 G=[-25 -4 0];


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
figure(1)
%axis equal
hold on
plot(0,0,'x','MarkerSize',50,'LineWidth',3)
visualizeOSM(Cam,osm_s);


%% filter - retrieve ojects in the field of view, in order !
%% compute projection
reprojection=computeProjection(Cam,osm_s);
visualizeProjection(reprojection,osm_s,Cam)


%% compute label area of projected images !
[labelLayers]=visualizeProjection_fillingLabels(reprojection,osm,Cam);
[prior]=mergeLayers(labelLayers);
figure
imagesc(prior)

outputfilename=strrep(outputfilename,'.xml','.mat');
outputfilename=strrep(outputfilename,'rawOSM','priors');
save(outputfilename,'prior')



% I=imread('0000000010.png');
% figure
% subplot(311)
% imagesc(prior);
% subplot(312)
% imshow(I);
% subplot(313)
% rgb=visualize_combinePriorIm(prior,I);
% imshow(rgb);
