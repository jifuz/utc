function main

%%%%%%%%%%%%%%%%%%%%
% input:
% -lat/lon
% -camera param (including heading)
%
% output:
% -map parsed from OpenStreetMap
%%%%%%%%%%%%%%%%%%%%

global sub_pathName
global GPS
global yawData
global fnames_test
global fnames
global osm

%% directory to dataset
sub_pathName='2011_09_30_drive_0028';
path_to_data=['data/' sub_pathName '/oxts/data'];


%% read data
path_to_names=['data/' sub_pathName '/images/test.txt'];
fnames_test=dlmread(path_to_names);

datafiles=dir([path_to_data '/*.txt']);
numfiles=length(datafiles);


%% read GPS
xs=[];
ys=[];
yawData=[];
fnames={};
for frame=0:numfiles-1
    
    if frame<1000
        n = sprintf('0000000%03d',frame);
    else
        n = sprintf('000000%03d',frame);
    end
    fnames{end+1}=n;
    
    if frame<1000
        fname = sprintf('%s/0000000%03d.txt',path_to_data,frame);
    else
        fname = sprintf('%s/000000%03d.txt',path_to_data,frame);
    end
    
    oxts_data = dlmread(fname);
    lat=oxts_data(1);
    lon=oxts_data(2);
    yawData(end+1)=oxts_data(6);

    xs(end+1)=lat;
    ys(end+1)=lon;
end
GPS=[xs' ys'];

%% estimate camera gaze direction
[GPS1_x,GPS1_y]=deg2utm(GPS(1,1),GPS(1,2));
[GPS2_x,GPS2_y]=deg2utm(GPS(2,1),GPS(2,2));
G=[GPS2_x-GPS1_x,GPS2_y-GPS1_y,0];


%% configure camera parameter
Cam=cameraParam(GPS(1,1),GPS(1,2),G); %Lat, Lon, Gaze direction


%% get map data from osm and parsing
outputfilename=sprintf(['data/' sub_pathName '/osm/OSM_%s_small.mat'],sub_pathName);
if ~exist(outputfilename,'file')
    padding=0.0002; %size of bounding box
    minlat = min(GPS(:,1))-padding;
    minlon = min(GPS(:,2))-padding;
    maxlat = max(GPS(:,1))+padding;
    maxlon = max(GPS(:,2))+padding;
    bboxIn=[minlon,minlat,maxlon,maxlat];
    outputfilename=sprintf(['data/' sub_pathName '/osm/OSM_%s_small.xml'],sub_pathName);
    
    %query osm and parsing
    osm=get_osm(bboxIn,outputfilename);
    
    if length(osm)==1 && isfield(osm,'osm')
        osm=osm.osm;
    end
    %save output (parsed)
    outputfilename=sprintf(['data/' sub_pathName '/osm/OSM_%s_small'],sub_pathName);
    save(outputfilename,'osm');
else
    osm=0;
    osm=load(outputfilename);
    osm=osm.osm;
end


%% specify visiblility level of each map object
global lowestLevel
global VisLevel
lowestLevel=10;
VisLevel=1;
[osm_s]=osmInfoStruct(osm);


%% visualize osm map data
figure(1)
%axis equal
hold on
plot(0,0,'x','MarkerSize',50,'LineWidth',3)
visualizeOSM(Cam,osm_s);


%% plot GPS trace
GPS_utm=GPS;
Gs=zeros(length(GPS),3); %sequence of camera gaze direction
for i=1:length(GPS)
    [lat lon]=deg2utm(GPS(i,1),GPS(i,2));
    GPS_utm(i,:)=[lat lon];
    Gs(i,1:2)=[cos(yawData(i)),sin(yawData(i))];
end
GPS_utm(2:end,1)=GPS_utm(2:end,1)-GPS_utm(1,1);
GPS_utm(2:end,2)=GPS_utm(2:end,2)-GPS_utm(1,2);
GPS_utm(1,:)=[0 0];

plot(GPS_utm(:,1),GPS_utm(:,2),'cx');

end