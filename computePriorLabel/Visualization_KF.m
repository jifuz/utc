

% figure
% ks=keys(reprojection)
% 
% for i=1:length(ks)
%     reprojection
% end

%% visualize complete data
figure
%axis equal
hold on
% plot(0,0,'x','MarkerSize',50,'LineWidth',3)
visualizeOSM(Cam,osm_s);


numfiles=340;
%% read timestamps
timestacmps=[];
path_to_time='../data/oxts/';
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
path_to_data='../data/oxts/data';
path_to_outputMap='map';


%function getMaps(path_to_data,path_to_outputMap,numfiles,zoom)

%path_to_outputMap='map/';
xs=[0];
ys=[0];
for frame=0:numfiles-1
    
    fname = sprintf('%s/0000000%03d.txt',path_to_data,frame);
%     outputFilename = sprintf('%s//0000000%03d.png',path_to_outputMap,frame);
    oxts_data = dlmread(fname);
    lat=oxts_data(1);
    lon=oxts_data(2);
    yawData=oxts_data(6);
    timestacmps(end+1)=oxts_data(6);
    locString=sprintf('%5.15f,%5.15f',lat,lon);
    
    
    [x,y]=deg2utm(lat, lon);
    if frame>0
        x=x-Cam.latUTM;
        y=y-Cam.lonUTM;
    else
        x=0;
        y=0;
    end
    
        fprintf('%d / %d \n',frame,numfiles+1);
        
        l=70;
        yyaw=-sin(yawData)*l;
        xyaw=cos(yawData)*l;
        hold on
        %plot([x xyaw],[y yyaw],'lineWidth',2);
        % pause(0.5)
    
        plot(x,y,'mx','LineWidth',2);
        
    xs(end+1)=x;
    ys(end+1)=y;
end
GPS=[xs' ys'];


%% load all ways
roads={};

nodes=osm_s.nodes;
ways=osm_s.ways;
nodesLevel=osm_s.nodesLevel;
waysLevel=osm_s.waysLevel;

ks=keys(osm_s.ways);

for i=1:ways.Count
   if waysLevel(ks{i})>VisLevel; %%not draw if level lower than visLevel
        continue;
    end
    way=ways(ks{i});   
    nds=way.nodes;
    x=zeros(numel(nds),1);
    y=zeros(numel(nds),1);
    for j=1:numel(nds)
        node=nodes(nds{j});
        [lat,lon]=deg2utm(str2num(node.LAT), str2num(node.LON));
        %transform in camera's coordinate
        lat=lat-Cam.latUTM;
        lon=lon-Cam.lonUTM;
        x(j)=lat;
        y(j)=lon;
    end
    
    
    tagS='';
    
    %%draw way tags
    tags=way.tag;
    fns=fieldnames(tags);
    
    residential=false;
    driveWay=false;
    service=false;
    footway=false;
  
    
    for tt=1:numel(fns)
        field=fns{tt};
        value=way.tag.(field);
        t=[field '=' value];
        tagS=[tagS '|' t];
        
        if numel(fns)==1
            field
        end
       
       if strcmp(lower(field), 'highway') ...
                && ( (strcmp(lower(value), 'residential')) ...
            || (strcmp(lower(value), 'living_street'))... 
            || (strcmp(lower(value), 'secondary')))
            residential=true;
       end
        
       if strcmp(lower(field), 'highway') ...
                && (strcmp(lower(value), 'tertiary'))
            driveWay=true;
       end
       
       if strcmp(lower(field), 'highway') ...
                && (strcmp(lower(value), 'service'))
            service=true;
       end
       
     
       
       if strcmp(lower(field), 'highway') ...
                && (strcmp(lower(value), 'footway')...
                || strcmp(lower(value), 'steps')...
                || strcmp(lower(value), 'access_ramp'))
            footway=true;
       end
        
     
       
    end
    
    
    firstNode=nodes(nds{1});
    lastNode=nodes(nds{numel(nds)});
    if residential || driveWay || service || footway
        roads{end+1,1}=[x y];       
    end
    
   
    
    
    
    
end



%% search for nearest road

roadsMat=cell2mat(roads);
roadsMatInd=zeros(size(roadsMat,1),1);
ind=1;
for i=1:length(roads)
    n=length(roads{i});
    roadsMatInd(ind:ind+n-1)=i;
    ind=ind+n;
    xmean=mean(roads{i}(:,1));
    ymean=mean(roads{i}(:,2));
    text(xmean,ymean,num2str(i))
end

numGPS=size(GPS,1);
countNN=zeros(numGPS,1);
for i=1:numGPS
    D=pdist2(GPS(i,:),roadsMat);
    [~,I]=min(D);
    countNN(i)=roadsMatInd(I);
end
roadInds=unique(countNN);
c=histc(countNN,sort(roadInds));
[count,I]=sort(c,1,'descend');
numRoadsConsider=2;
NNRoadInd=[];
NNRoad=[];
NNLane=[];
NNLaneCell={};
w=3;
for i=1:numRoadsConsider
    NNRoadInd(end+1)=roadInds(I(i));
    roadXY=roadsMat(roadsMatInd==NNRoadInd(end),:);
    [xl,yl,xr,yr]=roadBoundary(roadXY(:,1),roadXY(:,2),w);

    NNRoad=[NNRoad; roadXY];
    NNLane=[NNLane; [xr yr]];
    NNLaneCell{end+1}=[xr yr];
end

% plot(NNRoad(:,1),NNRoad(:,2),'yo');
% plot(NNLane(:,1),NNLane(:,2),'ys');



%% construct measurements for KF
yaw=zeros(size(GPS(:,1)));
dyaw=zeros(size(GPS(:,1)));
dx=zeros(size(GPS(:,1)));
dy=zeros(size(GPS(:,1)));
for i=2:length(yaw)
    vx=(GPS(i,1)-GPS(i-1,1));
    vy=(GPS(i,2)-GPS(i-1,2));
    angle=atan2(vy,vx);
    if angle<0
        angle=pi*2-angle;
    end
    yaw(i)=angle;
    
    dt=(time_s(i)-time_s(i-1));
    dyaw(i)=(yaw(i)-yaw(i-1))/dt;
    dx(i)=vx/dt;
    dy(i)=vy/dt;
end
yaw(1)=yaw(2);
dyaw(2)=0;
dyaw(1)=0;

measurements.x=GPS(:,1);
measurements.y=GPS(:,2);
measurements.yaw=yaw;
measurements.dx=dx;
measurements.dy=dy;
measurements.dyaw=dyaw;
measurements.time=time_s;
model.road=NNRoad;
model.lane=NNLane;
model.laneCell=NNLaneCell;

deltaT=time_s;
for i=1:length(time_s)-1
    deltaT(i)=time_s(i+1)-time_s(i);
end
deltaT(end)=deltaT(end-1);
model.time=deltaT;



[KF_xy,KF_dxdy,KF_V]=KF(measurements, model);
plot(KF_xy(:,1),KF_xy(:,2),'-','LineWidth',2);

for i=5:length(KF_xy)
    M=KF_xy(i,:);
    V=KF_V{i,1};
    C=V(1:2,1:2);
    [h,bp]=plot_gaussian_ellipsoid(M, C);
    plot(bp(1,:),bp(2,:),'y--');
end



%% 



