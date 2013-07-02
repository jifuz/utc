function reprojection=computeProjection(Cam,osm)

global VisLevel
global cutPlane
global buildingH
global buildingH_Var

reprojection=containers.Map; %projection:  ID->geometry
geometry=struct;
geometry.birdView={};
geometry.cameraView={};

nodes=osm.nodes;
ways=osm.ways;

%% compute projected ways geomerty
ks=keys(osm.ways);
for i=1:ways.Count
    
    %% filter those with level lower than vis threshold
    if osm.waysLevel(ks{i})>VisLevel
        continue;
    end
    
    %% get way geometry
    way=ways(ks{i});
    
    nds=way.nodes;
    x=zeros(numel(nds),1);
    y=zeros(numel(nds),1);
    for j=1:numel(nds)
        node=nodes(nds{j});
        [lat,lon]=deg2utm(str2num(node.LAT), str2num(node.LON));
        
        %transform in camera's coordinate
        lat=lat-Cam.latUTM-1.4;
        lon=lon-Cam.lonUTM-1.4;
        
        x(j)=lat;
        y(j)=lon;
    end
    
    
    %% process way tags
    tags=way.tag;
    
    
    
      %%%%%%rotate every thing w.r.t cam pointing to north
    [x_n,y_n]=rotateToNorth(x,y,Cam);
    x=x_n';
    y=y_n';
    
    %% compute projected images
    %filter geometry outside of camera view!!!!
    if sum(y<cutPlane)>0
        %filter geometry outside of camera view, denote by ind(i)=1
%         [xVis,yVis,ind]=VisFilter(Cam,x,y);
%         
%         if (strcmp('highway',classifyType(tags)) || strcmp('footway',classifyType(tags)))
%             [xVis,yVis,ind]=VisFilter_road(Cam,x,y);
%         end
        
        [xVis,yVis,ind]=VisFilter_road(Cam,x,y);
        %filter those not in the scene
        if sum(ind)==length(ind)
            continue;
        end
    
         x=xVis;
         y=yVis;
    end
    
    %in case for road boundary
    if strcmp('highway',classifyType(tags)) || strcmp('footway',classifyType(tags))
        w=roadWidth(lower(tags.('highway')));
        [xl,yl,xr,yr]=roadBoundary(x,y,w);
        [xl,yl,ind]=VisFilter_road(Cam,xl,yl);
        [xr,yr,ind]=VisFilter_road(Cam,xr,yr);
        [xl_n,yl_n]=rotateToNorthInverse(xl,yl,Cam);
        [xr_n,yr_n]=rotateToNorthInverse(xr,yr,Cam);
        xl=xl_n';yl=yl_n';
        xr=xr_n';yr=yr_n';
    end
    
    [x_n,y_n]=rotateToNorthInverse(x,y,Cam);
    x=x_n';
    y=y_n';
%     %filter geometry outside of camera view, denote by ind(i)=1
%     [xVis,yVis,ind]=VisFilter(Cam,x,y);
%     
%     %filter those not in the scene
%     if sum(ind)==length(ind)
%         continue;
%     end
%     
%     %cut highway point out of camera view to avoid junk aritfacts
%     if (strcmp('highway',classifyType(tags)) || strcmp('footway',classifyType(tags))) && sum(ind)>0 %if there are lines across camera plane
%        
%         for ii=1:numel(ind)-1
%             if ind(ii)~=ind(ii+1)
%                 %compute the point intersect with camera screen just before
%                 %it goes behind the screen plane
%                 a=[x(ii) y(ii)];
%                 b=[x(ii+1) y(ii+1)];
%                 [xx,yy]=pointsOnCameraPlane(a,b,Cam);
%                 break;
%             end
%         end
%         
%         
%         if ind(ii)==1 %first half are 1, behind screen, cut
%             
%               x=[xx; xVis(ii+1:end)];
%               y=[yy; yVis(ii+1:end)];
% %                        x=xVis(ii+1:end);
% %                        y=yVis(ii+1:end);
%         else %last half are 1, behind screen, cut
%               x=[xVis(1:ii); xx];
%               y=[yVis(1:ii); yy];
% %                        x=xVis(1:ii);
% %                        y=yVis(1:ii);
%         end
%         
%         %non highway points don't need to be cut since artifacts are trivial
%     else
%         x=xVis;
%         y=yVis;
%     end
    
    
    
    
    %compute projected geometry in the camera view
    [xp,yp]=projection_KFSample(Cam,x,y,zeros(size(x,1),size(x,2)));
    
    %     %filter those not in the scene
    %     hasPointInTheScene=false;
    %     for ii=1:numel(xp)
    %         if xp(ii)>0 && xp(ii)<Cam.nx && yp(ii)>0 && yp(ii)<Cam.ny
    %             hasPointInTheScene=true;
    %         end
    %     end
    %     if ~hasPointInTheScene
    %         continue;
    %     end
    
    if strcmp('building',classifyType(tags)) %building bottom/roof
        h=mvnrnd(buildingH,buildingH_Var,1);
        [xptop,yptop]=projection_KFSample(Cam,x,y,h*ones(size(x,1),size(x,2)));
        geometry.cameraView{1}=[xp yp];
        geometry.cameraView{2}=[xptop yptop];
        
    elseif (strcmp('highway',classifyType(tags)) || strcmp('footway',classifyType(tags))) && length(x)>1%highway, left/right boundary
%         %retrieve width of highway class
%         w=roadWidth(lower(tags.('highway')));
%         %compute boundary
%         [xl,yl,xr,yr]=roadBoundary(x,y,w);
%         xl0=xl;
%         yl0=yl;
%         [xl,yl,ind]=VisFilter(Cam,xl,yl);
% %         sum(abs(xl-xl0))
% %         sum(abs(yl-yl0))
%         
%         [xr,yr,ind]=VisFilter(Cam,xr,yr);
        %compute reprojected geometry
        [xpl,ypl]=projection_KFSample(Cam, xl,yl,zeros(size(xl,1),size(xl,2)));
        [xpr,ypr]=projection_KFSample(Cam, xr,yr,zeros(size(xr,1),size(xr,2)));
        
        geometry.birdView{2}=[xl yl];
        geometry.birdView{3}=[xr yr];
        geometry.cameraView{1}=[xp yp];
        geometry.cameraView{2}=[xpl ypl];
        geometry.cameraView{3}=[xpr ypr];
        
    else  %other map element
        geometry.cameraView{1}=[xp yp];
    end
    
    %store the projected geometry
    geometry.birdView{1}=[x y];
    reprojection(ks{i})=geometry;
end





%% compute projected nodes geomerty
nodeBoxW=6;
ks=keys(osm.nodes);
for i=1:nodes.Count
  
%% filter those with level lower than vis threshold
    if osm.nodesLevel(ks{i})>VisLevel
        continue;
    end
    
    %% get node location, and compute bounding box
    node=nodes(ks{i});
    [lat,lon]=deg2utm(str2num(node.LAT), str2num(node.LON));
    xc=lat-Cam.latUTM;
    yc=lon-Cam.lonUTM;
    w=nodeBoxW/2;
    x=[xc+w; xc-w; xc-w; xc+w; xc+w;];
    y=[yc+w; yc+w; yc-w; yc-w; yc+w;];
    
    
    %% process way tags
    if ~isfield(node,'tag')
        continue;
    end
    
    tags=node.tag;

    %% compute projected images    
    %filter geometry outside of camera view, denote by ind(i)=1
    [xVis,yVis,ind]=VisFilter(Cam,x,y);
    
    %filter those not in the scene
    if sum(ind)==length(ind)
        continue;
    end
    
    %compute projected geometry in the camera view
%     tic
%     [xp,yp]=projection_KFSample(Cam,x,y,zeros(size(x,1),size(x,2)));
%     tt=toc
%     tic
    [xp,yp]=projection_KFSample(Cam,x,y,zeros(size(x,1),size(x,2)));
%     ttt=toc
%     ttt/tt
    %store the projected geometry
    geometry.cameraView{1}=[xp yp];
    geometry.birdView{1}=[x y];
    reprojection(ks{i})=geometry;
    
end



end



function w=roadWidth(type)

switch lower(type)
    case 'residential'
        w=6;
    case 'tertiary'
        w=15;
    case 'service'
        w=10;%w=4;
    case 'secondary'
        w=4;%w=4;
    case 'footway'
        w=2;%w=4;
    otherwise %maybe mostly footway?
        w=10;%w=2;
end

end



% intersection of line a-b and line camera vertG
function [xx,yy]=pointsOnCameraPlane(a,b,Cam)

% vertG=vertVector(Cam.G);
%
% %y=a0*x, line of vertG
% b0=0;
% a0=vertG(2)/vertG(1);
%
% %y=a1*x+b1, line of a and b
% a1=(b(2)-a(2))/(b(1)-a(1));
% b1=b(2)-a1*b(1);
%
%
% xx=-(b0-b1)/(a0-a1);
% yy=a0*xx+b0;

[ar,br]=rotateToNorth(a,b,Cam);

% yOnPlane=-Cam.n+0.1;
 
yOnPlane=-Cam.n;
r=(yOnPlane-ar(2))/(br(2)-ar(2));
% xOnPlane=ar(1)+r*(br(1)-ar(1));
% 
% [xx,yy]=rotateToNorthInverse(xOnPlane,yOnPlane,Cam);

xx=a(1)+r*(b(1)-a(1));
yy=a(2)+r*(b(2)-a(2));

end

function vertG=vertVector(G)
vertG=-[G(2) -G(1)];
end