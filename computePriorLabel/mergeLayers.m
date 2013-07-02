function [prior]=mergeLayers(labelLayers)

%label mapping, 1 is at bottom layer
label=containers.Map;
label('1')='parking';%gravel floor
label('2')='footway';%gravel floor
label('3')='highway';%asphalt floor
label('4')='building';%building
label('5')='transportation_objects';%object

ks=keys(labelLayers);
prior=zeros(size(labelLayers(ks{1})));

for i=1:label.Count
    if ~isKey(labelLayers,label(num2str(i)))
        continue;
    end
    currLayer=labelLayers(label(num2str(i)));
    prior(currLayer==1)=i;
end

% prior(prior==1)=gravel_floor_ind;
% prior(prior==2)=gravel_floor_ind;
% prior(prior==3)=asphalt_floor_ind;
% prior(prior==4)=building_ind;
% prior(prior==5)=object_ind;




end