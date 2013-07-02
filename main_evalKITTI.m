global KF_Qnoise
global KF_Rnoise
global cutPlane
KF_Qnoise=3;
KF_Rnoise=3;
cutPlane=KF_Qnoise*2;%visFilter cutoff plane
global buildingH
global buildingH_Var
buildingH=60;
buildingH_Var=60;
global numSample

close all

%% main
dname='2011_09_26_drive_0059';
dname='2011_10_03_drive_0034';
sub_pathName='2011_09_30_drive_0028_image_03';
sub_pathName='2011_10_03_drive_0034_image_03';
path_to_time=['data/' sub_pathName '/oxts/'];
path_to_data=['data/' sub_pathName '/oxts/data'];


%% read data

path_to_names=['data/' sub_pathName '/images/test.txt'];
fnames_test=dlmread(path_to_names);

datafiles=dir([path_to_data '/*.txt']);
numfiles=length(datafiles);
%% read timestamps
timestacmps=[];
timefname=[path_to_time 'timestamps.txt'];
time_data = dlmread(timefname);
time_s=zeros(length(time_data),1);
t0=time_data(1,:);
for i=2:length(time_data)
    line=time_data(i,:);
    time_s(i)=(line(2)-t0(2))*60+(line(3)-t0(3));
end
time_s=[time_s; time_s(end)+(time_s(end)-time_s(end-1))];

%% read GPS
%path_to_outputMap='map';

%function getMaps(path_to_data,path_to_outputMap,numfiles,zoom)

%path_to_outputMap='map/';
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
    
%     outputFilename = sprintf('%s//0000000%03d.png',path_to_outputMap,frame);
    oxts_data = dlmread(fname);
    lat=oxts_data(1);
    lon=oxts_data(2);
    yawData(end+1)=oxts_data(6);
    timestacmps(end+1)=oxts_data(6);
    locString=sprintf('%5.15f,%5.15f',lat,lon);
    
    xs(end+1)=lat;
    ys(end+1)=lon;
end
GPS=[xs' ys'];

[GPS1_x,GPS1_y]=deg2utm(GPS(1,1),GPS(1,2));
[GPS2_x,GPS2_y]=deg2utm(GPS(2,1),GPS(2,2));
G=[GPS2_x-GPS1_x,GPS2_y-GPS1_y,0];
Cam=cameraParam(GPS(1,1),GPS(1,2),G);


%% get all osm
outputfilename=sprintf(['data/' sub_pathName '/osm/OSM_%s_small.mat'],dname);
if exist(outputfilename,'file')==0
    padding=0.0002;
    minlat = min(GPS(:,1))-padding;
    minlon = min(GPS(:,2))-padding;
    maxlat = max(GPS(:,1))+padding;
    maxlon = max(GPS(:,2))+padding;
    bboxIn=[minlon,minlat,maxlon,maxlat];
    outputfilename=sprintf(['data/' sub_pathName '/osm/OSM_%s_small.xml'],dname);
    %query osm
    osm=get_osm(bboxIn,outputfilename);
    if length(osm)==1 && isfield(osm,'osm')
        osm=osm.osm;
    end
    %save output (parsed)
    outputfilename=sprintf(['data/' sub_pathName '/osm/OSM_%s_small'],dname);
%     save(outputfilename,'osm');
else
    osm=0;
    osm=load(outputfilename);
    osm=osm.osm;
end

% specify information structure - visiblility level of each obj
global lowestLevel
global VisLevel
lowestLevel=10;
VisLevel=1;
[osm_s]=osmInfoStruct(osm);


% visualize complete data
figure(1)
%axis equal
hold on
plot(0,0,'x','MarkerSize',50,'LineWidth',3)
visualizeOSM(Cam,osm_s);

GPS_utm=GPS;
Gs=zeros(length(GPS),3);
for i=1:length(GPS)
    [lat lon]=deg2utm(GPS(i,1),GPS(i,2));
    GPS_utm(i,:)=[lat lon];
    Gs(i,1:2)=[cos(yawData(i)),sin(yawData(i))]; 
end
GPS_utm(2:end,1)=GPS_utm(2:end,1)-GPS_utm(1,1);
GPS_utm(2:end,2)=GPS_utm(2:end,2)-GPS_utm(1,2);
GPS_utm(1,:)=[0 0];

%      [x_n,y_n]=rotateToNorth(GPS_utm(:,1),GPS_utm(:,2),Cam);
%      GPS_utm(:,1)=x_n';
%      GPS_utm(:,2)=y_n';

plot(GPS_utm(:,1),GPS_utm(:,2),'cx');
for i=1:length(GPS)
    currX_v=[GPS_utm(i,1); GPS_utm(i,1)+Gs(i,1)*0.4];
    currY_v=[GPS_utm(i,2); GPS_utm(i,2)+Gs(i,2)*0.4];
%     [x_n,y_n]=rotateToNorth(currX_v,currY_v,Cam);
%     currX_v=x_n';
%     currY_v=y_n';
%     plot(currX_v,currY_v,'-');
end


%% KF all
[KF_xy,KF_dxdy,KF_V]=KF_main(GPS,time_s,osm_s);
% [KF_xy,KF_dxdy,KF_V]=KF(measurements, model);
% plot(KF_xy(:,1),KF_xy(:,2),'-','LineWidth',2);
% 
% for i=5:length(KF_xy)
%     M=KF_xy(i,:);
%     V=KF_V{i,1};
%     C=V(1:2,1:2);
%     [h,bp]=plot_gaussian_ellipsoid(M, C);
%     plot(bp(1,:),bp(2,:),'y--');
% end
C=KF_V;
C_dxdy=C;
for i=1:length(KF_V)
     VV=KF_V{i,1};
     C{i,1}=VV(1:2,1:2);
     C_dxdy{i,1}=VV(3:4,3:4);
%     [h,bp]=plot_gaussian_ellipsoid(M, C);
%     plot(bp(1,:),bp(2,:),'y--');
% end
end

numSample=4;
global KFpert
global currKFpert
KFpert=cell(numSample,1);
for i=1:length(KF_xy)
    M=KF_xy(i,:);
    xy=mvnrnd(M,C{i,1},numSample);
    dxdy=mvnrnd(KF_dxdy(i,:),C_dxdy{i,1},numSample);
    dxdy=dxdy./norm(dxdy);
%     xy(:,1)-M(1);
    KFpert{i}=[xy(:,1)-M(1) xy(:,2)-M(2) dxdy(:,1) dxdy(:,2)];
end

%% compute priors
priors.fnames={};
priors.priors={};
for ii=29:length(fnames_test)
    
    currTest=fnames_test(ii);
    currKFpert=KFpert{currTest};
    i=currTest+1;
    
    if ii>1 %update cam parameter
        [GPS1_x,GPS1_y]=deg2utm(GPS(i-1,1),GPS(i-1,2));
        [GPS2_x,GPS2_y]=deg2utm(GPS(i,1),GPS(i,2));
        G=[GPS2_x-GPS1_x,GPS2_y-GPS1_y,0];
        Cam=cameraParam(GPS(i,1),GPS(i,2),Gs(i,:));
    end
    
    fname=fnames{i};
    echo on
        fprintf('processing prior %s (%d/%d)...\n',fname,ii,length(fnames_test));
echo off
    % compute projection
    tic
  reprojections=computeProjection_KFSample(Cam,osm_s);
     reprojection=computeProjection(Cam,osm_s);

    tt=toc
    visualizeProjection(reprojection,osm_s,Cam)
    
    
    % compute label area of projected images !
%     [labelLayers]=visualizeProjection_fillingLabels_KFSample(reprojections,osm,Cam);
 %    [labelLayers]=visualizeProjection_fillingLabels(reprojection,osm,Cam);
       [labelLayers,priorLayers]=visualizeProjection_fillingLabels_KFSample(reprojections,osm,Cam);

    %[prior]=mergeLayers(labelLayers);
    

%     prior=blurPrior(C{ii},prior);
%     [v,d]=eig(C); 

    
    %visualize
    %     figure
    %     imagesc(prior);axis image;
    
    %save prior jpg and mat to disk
    prior_outputfilename=sprintf(['data/' sub_pathName '/priors/%s'],fname);
    saveas(gcf,prior_outputfilename,'jpg');
    save(prior_outputfilename,'priorLayers','labelLayers')
    priors.fnames{end+1}=fname;
    priors.priors{end+1}=priorLayers;
    
    %save the overlaid prior/image for visualization
    imName=sprintf(['data/' sub_pathName '/images/%s.jpg'],fname);
    im=imread(imName);
    im_overlaid=overlayPrior(labelLayers,im);
    figure
    imagesc(im_overlaid)
    imOverlaidName=sprintf(['data/' sub_pathName '/priors/%s_overlaid.jpg'],fname);
    axis equal
%     saveas(gcf,imOverlaidName,'jpg');
    imwrite(im_overlaid,imOverlaidName,'jpg')
     
    close all
    fprintf('DONE processing prior %s (%d/%d)...\n',fname,ii,length(fnames_test));
    


% prior_outputfilename=sprintf('data/sample/priors/%s',fname);
% 
%         priorLayers=load(prior_outputfilename);
% %save the overlaid prior/image for visualization
%     imName=sprintf('data/sample/images/%s.jpg',fname);
%     im=imread(imName);
%     im_overlaid=overlayPrior(labelLayers,im);
%     figure
%     imagesc(im_overlaid)
%     imOverlaidName=sprintf('data/sample/priors/%s_overlaid.jpg',fname);
%     axis equal
% %     saveas(gcf,imOverlaidName,'jpg');
%     imwrite(im_overlaid,imOverlaidName,'jpg')
%          close all

end
save(['data/' sub_pathName '/priors/priors'],'priors');


return

%% combine - maxent 
dataDir='data/sample/images/';
inferDir='data/sample/images/ex_infer/';
groundTruthDir='data/sample/images/labels/';
priorDir='data/sample/priors';
listFile='test.txt';
alpha=0.5;


dataDir='./test_data/';
inferDir='./test_data/ex_infer/';
groundTruthDir='./test_data/labels/';
priorDir='./test_data/priors';
listFile='test.txt';

infer_test_func(dataDir,inferDir,groundTruthDir,priorDir,listFile,alpha)




%% evaluation


%% 
% fname='00000000xx';
% im=imread([imDir '/' fname '.jpg']);
% 





%% 


