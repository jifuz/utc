function [prior]=mergeLayers_Blur(labelLayers)

%label mapping, 1 is at bottom layer
label=containers.Map;
label('1')='parking';
label('2')='footway';
label('3')='highway';
label('4')='building';
label('5')='transportation_objects';

ks=keys(labelLayers);
prior=zeros(size(labelLayers(ks{1})));

for i=1:label.Count
    if ~isKey(labelLayers,label(num2str(i)))
        continue;
    end
    currLayer=labelLayers(label(num2str(i)));
    
    blurCurrLayer=blur(currLayer);
    
    prior(currLayer==1)=i;
end
end

function blurCurrLayer=blur(currLayer)

blurredImage = conv2(redChannel, ones(windowWidth)/windowWidth^2, 'same');

gaussianFilter = fspecial('gaussian', [7, 7], 5);
gaussianPomegranate = imfilter(pomegranate, gaussianFilter, 'symmetric', 'conv');
subplot(1,2,1), image(pomegranate);
subplot(1,2,2), image(gaussianPomegranate), title('Blurred Pomegranate, blur matrix size 7');

end