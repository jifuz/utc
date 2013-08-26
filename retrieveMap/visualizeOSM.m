function visualizeOSM(Cam,osm)

global VisLevel


% [minlat minlon]=deg2utm(str2num(osm.bounds.minlat),str2num(osm.bounds.minlon));
[maxlat maxlon]=deg2utm(str2num(osm.bounds.maxlat),str2num(osm.bounds.maxlon));



nodes=osm.nodes;
ways=osm.ways;
nodesLevel=osm.nodesLevel;
waysLevel=osm.waysLevel;

%% plot nodes
ks=keys(osm.nodes);
for i=1:nodes.Count
    if nodesLevel(ks{i})>VisLevel; %%not draw if level lower than visLevel
        continue;
    end
    
    node=nodes(ks{i});
    lat=str2num(node.LAT);
    lon=str2num(node.LON);
    [lat,lon]=deg2utm(lat,lon);
    %transform in camera's coordinate
    lat=lat-Cam.latUTM;
    lon=lon-Cam.lonUTM;
    if isfield(node,'tag')
        fns=fieldnames(node.tag);
        field=fns{1};
        value=node.tag.(field);
        if strcmp(lower(field), 'barrier')...
        || strcmp(lower(field), 'amenity')...        
        || strcmp(lower(field), 'entrance')        

        else
        %text(lat,lon,[field '=' value]);
        end
        plot(lat,lon,'o')
    end
    
end


%% plot ways
ks=keys(osm.ways);

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
    
    %%%%%%rotate every thing w.r.t cam pointing to north
%     [x_n,y_n]=rotateToNorth(x,y,Cam);
%     x=x_n';
%     y=y_n';
     
    tagS='';
    
    %%draw way tags
    tags=way.tag;
    fns=fieldnames(tags);
    
    building=false;
    residential=false;
    driveWay=false;
    service=false;
    parking=false;
    footway=false;
    nolabel=false;
    barrier=false;
    park=false;
    
    for tt=1:numel(fns)
        field=fns{tt};
        value=way.tag.(field);
        t=[field '=' value];
        tagS=[tagS '|' t];
        
        if numel(fns)==1
            field
        end
        if strcmp(lower(field), 'building') %...
%                 && (strcmp(lower(value), 'yes')...
%                 || strcmp(lower(value), 'true'))
            building=true;
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
       
       if strcmp(lower(field), 'amenity') ...
                && (strcmp(lower(value), 'parking'))
            parking=true;
       end
       
        
       if (strcmp(lower(field), 'barrier'))
            barrier=true;
       end
       
       if (strcmp(lower(field), 'amenity') && (strcmp(lower(value), 'post_box')))
            nolabel=true;
       end
       
       if strcmp(lower(field), 'highway') ...
                && (strcmp(lower(value), 'footway')...
                || strcmp(lower(value), 'steps')...
                || strcmp(lower(value), 'access_ramp'))
            footway=true;
       end
        
         if strcmp(classifyType(way.tag),'park')
            park=true;
         end
       
    end
    
    
    firstNode=nodes(nds{1});
    lastNode=nodes(nds{numel(nds)});
    if building %&&strcmp(lastNode.ID, firstNode.ID)
        
%fill(x,y,'r');
        drawPolyline([x,y],'lineWidth', 6,'Color','r')
%     elseif strcmp(lastNode.ID, firstNode.ID) && (~building) %%area
%         drawPolyline([x,y],'lineWidth', 1,'LineStyle','y')
    elseif residential
        drawPolyline([x,y],'lineWidth', 4,'Color','k');

    elseif driveWay
        drawPolyline([x,y],'lineWidth', 10,'Color','k');
    elseif service
        drawPolyline([x,y],'lineWidth', 2,'Color','k');
    elseif parking
        drawPolyline([x,y],'lineWidth', 2,'Color','y');
    elseif footway
        drawPolyline([x,y],'lineStyle', '--','Color','k');
    elseif barrier
        drawPolyline([x,y],'lineStyle', 'x');
        text(mean(x),mean(y),[tagS '..barrier']);
    elseif park
        drawPolyline([x,y],'lineWidth', 2,'Color','g');
        
    elseif nolabel
        drawPolyline([x,y]);
    else
        drawPolyline(x,y);
        text(mean(x),mean(y),[tagS '...else:' ks{i} '|' num2str(numel(fns))]);

    end
    
   
    
end