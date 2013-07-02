function osm=get_osm(bboxIn,outputfilename)


%osm=urlread('http://api.openstreetmap.org/api/0.6/map?bbox=11.54,48.14,11.543,48.145');

%Forbes Smith
%40.444767,-79.946995
%40.444089,-79.944838

%Forbes-Murray
40.439133,-79.925097
40.437328,-79.921331

minLon=bboxIn(1);
minLat=bboxIn(2);
maxLon=bboxIn(3);
maxLat=bboxIn(4);


bbox=sprintf('%s,%s,%s,%s',num2str(minLon,10),num2str(minLat,10),num2str(maxLon,10),num2str(maxLat,10));
%url='http://api.openstreetmap.org/api/0.6/map?bbox=-79.946995,40.444089,-79.944838,40.444767';
%url='http://api.openstreetmap.org/api/0.6/map?bbox=-79.925097,40.437328,-79.921331,40.439133';

base='http://api.openstreetmap.org/api/0.6/map?bbox=';
url=[base bbox];

%url='http://api.openstreetmap.org/api/0.6/map?bbox=11.54,48.14,11.541,48.141';
%filename='Forbes_Murray1.xml';
filename=outputfilename;
urlwrite(url,filename);

%osm_DOM=xmlread(filename);
osm = xmltools(filename);
osm=parseOSM(osm);
end