function [osm_s]=osmInfoStruct(osm)

%%
% specify what detail level the map should display, the higher the more
% visible
%%

global lowestLevel
spec=infoSpec(); % spec the level

%osm_f=osm;
osm_s=osm;


nodes=osm.nodes;
ways=osm.ways;
wlevel=containers.Map;
nlevel=containers.Map;


%% ways level
ks=keys(osm.ways);
for i=1:ways.Count
    l=lowestLevel;
    way=ways(ks{i});
    tags=way.tag;
    fns=fieldnames(tags);
    for tt=1:numel(fns)
        field=fns{tt};
        value=way.tag.(field);
        if spec.isKey(lower(field))
            level=spec(lower(field));
            if isnumeric(level)
                l=level;
            elseif level.isKey(lower(value))
                l=level(lower(value));
            end
        end
        
        
    end
    %%assign way level
    wlevel(ks{i})=l;
    %%assign element node to the same level
    wayNodes=way.nodes;
    for nn=1:numel(wayNodes)
        nlevel(wayNodes{nn})=l;
    end
    
%     if l>VisLevel %filter out objects of vislevel<threshold
%         remove(osm_f.ways,ks{i});
%     end
end



%% nodes level
ks=keys(osm.nodes);
for i=1:nodes.Count
    l=lowestLevel;
    node=nodes(ks{i});
    
    %handle node that has no tag 
    if ~isfield(node,'tag') 
        if ~nlevel.isKey(ks{i}) %if not an element of a way,lowest level
            nlevel(ks{i})=l;
        end
        continue;
    end
    
    %handle node that has tag 
    tags=node.tag;
    fns=fieldnames(tags);
    for tt=1:numel(fns)
        field=fns{tt};
        value=node.tag.(field);
        if spec.isKey(lower(field))
            level=spec(lower(field));
            if isnumeric(level)
                l=level;
            elseif level.isKey(lower(value))
                l=level(lower(value));
            end
        end
    end
    
    %%assign level 
    %%if it is not part of a way (in which case = way-node level)
    %%or if node level is higher than assigned way-node level
    if ~nlevel.isKey(ks{i}) || (l<nlevel(ks{i}))
        nlevel(ks{i})=l;
    end
    
%     if l>VisLevel %filter out objects of vislevel<threshold
%         remove(osm_f.nodes,ks{i});
%     end
end


osm_s.nodesLevel=nlevel;
osm_s.waysLevel=wlevel;

end

%% specify the level, the higher (1) the more visible
function spec=infoSpec()

spec=containers.Map;

spec('building')=1;
spec('barrier')=3;

amenityLevel=containers.Map;
amenityLevel('parking')=2;
amenityLevel('post_box')=3;
spec('amenity')=amenityLevel;

highwayLevel=containers.Map;
highwayLevel('residential')=1;
highwayLevel('living_street')=1;
highwayLevel('secondary')=1;

highwayLevel('tertiary')=1;
highwayLevel('service')=1;
highwayLevel('bus_stop')=1;
highwayLevel('traffic_signals')=1;
highwayLevel('footway')=1;
highwayLevel('steps')=1;
highwayLevel('access_ramp')=1;

spec('highway')=highwayLevel;


landuseLevel=containers.Map;
landuseLevel('recreation_ground')=1;
spec('landuse')=landuseLevel;

leisureLevel=containers.Map;
leisureLevel('park')=3;
spec('leisure')=leisureLevel;


end
