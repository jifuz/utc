function [osm,outputfilename]=retrieveData(lat, lon,G)

latDiff=0.0005;
lonDiff=0.002;
c=clock;

%outputfilename=sprintf('%s-%s-%s_tempOSM.xml',num2str(c(1)),num2str(c(2)),num2str(c(3)));

minlon=lon-lonDiff;
minlat=lat-latDiff;
maxlon=lon+lonDiff;
maxlat=lat+latDiff;
bboxIn=[minlon,minlat,maxlon,maxlat];

outputfilename=sprintf('../data/rawOSM/OSM_%s_%s__%s-%s-%s.xml',...
    num2str(lat),num2str(lon),...
    num2str(c(1)),num2str(c(2)),num2str(c(3)));

% parsedOSMFileName=outputfilename;
% parsedOSMFileName=strrep(parsedOSMFileName,'.xml','.mat');
% parsedOSMFileName=strrep(parsedOSMFileName,'rawOSM','parsedOSM');
% if exist(parsedOSMFileName, 'file')
%     osm=load(parsedOSMFileName);
%     osm=osm.osm;
% else    
    osm=get_osm(bboxIn,outputfilename);
%     save(parsedOSMFileName,'osm');
% end
    
end