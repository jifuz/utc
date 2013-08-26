function osm=get_osm(bboxIn,outputfilename)

%%
% retrive map of bounding box - bboxIn
%%

%osm=urlread('http://api.openstreetmap.org/api/0.6/map?bbox=11.54,48.14,11.543,48.145');

minLon=bboxIn(1);
minLat=bboxIn(2);
maxLon=bboxIn(3);
maxLat=bboxIn(4);

bbox=sprintf('%s,%s,%s,%s',num2str(minLon,10),num2str(minLat,10),num2str(maxLon,10),num2str(maxLat,10));
base='http://api.openstreetmap.org/api/0.6/map?bbox=';
url=[base bbox];

filename=outputfilename;
urlwrite(url,filename);

%osm_DOM=xmlread(filename);
osm = xmltools(filename); %xml to matlab
osm=parseOSM(osm); %parse to "hashmap" structure
end