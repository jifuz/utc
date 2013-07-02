function reprojection=computeProjection(Cam,osm)

global VisLevel

buildingH=160;

numSample=50;
reprojections=cell(numSample);
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
    
    
    %% compute projected images
    
    %filter geometry outside of camera view, denote by ind(i)=1
    [xVis,yVis,ind]=VisFilter(Cam,x,y);
    
    %filter those not in the scene
    if sum(ind)==length(ind)
        continue;
    end
    
    %cut highway point out of camera view to avoid junk aritfacts
    if (strcmp('highway',classifyType(tags)) || strcmp('footway',classifyType(tags))) && sum(ind)>0 %if there are lines across camera plane
        
        for ii=1:numel(ind)-1
            if ind(ii)~=ind(ii+1)
                %compute the point intersect with camera screen just before
                %it goes behind the screen plane
                a=[x(ii) y(ii)];
                b=[x(ii+1) y(ii+1)];
                [xx,yy]=pointsOnCameraPlane(a,b,Cam);
                break;
            end
        end
        
        
        if ind(ii)==1 %first half are 1, behind screen, cut
            
            x=[xx; xVis(ii+1:end)];
            y=[yy; yVis(ii+1:end)];
            %                        x=xVis(ii+1:end);
            %                        y=yVis(ii+1:end);
        else %last half are 1, behind screen, cut
            x=[xVis(1:ii); xx];
            y=[yVis(1:ii); yy];
            %                        x=xVis(1:ii);
            %                        y=yVis(1:ii);
        end
        
        %non highway points don't need to be cut since artifacts are trivial
    else
        x=xVis;
        y=yVis;
    end
    
    
    
    
    %compute projected geometry in the camera view
    tic
    [xp,yp]=projection(Cam,x,y,zeros(size(x,1),size(x,2)));
    
    for iii=1:numSample
        x_=x+rand(1);
        y_=y+rand(1);
        [xp,yp]=projection(Cam,x_,y_,zeros(size(x,1),size(x,2)));
        %     end
        %     toc
        
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
            [xptop,yptop]=projection(Cam,x,y,buildingH*ones(size(x,1),size(x,2)));
            geometry.cameraView{1}=[xp yp];
            geometry.cameraView{2}=[xptop yptop];
            
        elseif (strcmp('highway',classifyType(tags)) || strcmp('footway',classifyType(tags))) && length(x)>1%highway, left/right boundary
            %retrieve width of highway class
            w=roadWidth(lower(tags.('highway')));
            %compute boundary
            [xl,yl,xr,yr]=roadBoundary(x,y,w);
            xl0=xl;
            yl0=yl;
            [xl,yl,ind]=VisFilter(Cam,xl,yl);
%             sum(abs(xl-xl0));
%             sum(abs(yl-yl0));
%             
            [xr,yr,ind]=VisFilter(Cam,xr,yr);
            %compute reprojected geometry
            [xpl,ypl]=projection(Cam, xl,yl,zeros(size(x,1),size(x,2)));
            [xpr,ypr]=projection(Cam, xr,yr,zeros(size(x,1),size(x,2)));
            
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
        reprojections{iii}=reprojection;%%KF
    end
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
    [xp,yp]=projection(Cam,x,y,zeros(size(x,1),size(x,2)));
    
    for iii=1:numSample
        x_=x+rand(1);
        y_=y+rand(1);
        [xp,yp]=projection(Cam,x_,y_,zeros(size(x,1),size(x,2)));
        
        %store the projected geometry
        geometry.cameraView{1}=[xp yp];
        geometry.birdView{1}=[x y];
        reprojection=reprojections{iii};%%KF
        reprojection(ks{i})=geometry;
        reprojections{iii}=reprojection;%%KF
    end
end



end



function w=roadWidth(type)

switch lower(type)
    case 'residential'
        w=12;
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