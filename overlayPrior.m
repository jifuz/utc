function rgb=overlayPrior(labelLayers,sourceImage)




%make sure building occludes roads!!!
layers={'building','highway'};
% buildingLayer=labelLayers('building');
% highwayLayer=labelLayers('highway');
% highwayLayer(buildingLayer>0.8)=0;
% labelLayers('highway')=highwayLayer;

layerColor=containers.Map;
layerColor('highway')=[256 256 0]; %asphalt
layerColor('building')=[256 0 0];
layerColor('footway')=[256 0 0]; %concrete
layerColor('parking')=[256 0 0]; %gravel
layerColor('park')=[256 0 0];
layerPrior=containers.Map;

colorPrior=zeros(240,795,3);
for i=1:length(layers)
    currLayer=layers{i};
    ins=labelLayers(currLayer);
    ins=ins';
    %ins(ins<0.4)=0;
    layerPrior(currLayer)=ins;
    color=layerColor(currLayer);
    R=color(1).*ins;
    G=color(2).*ins;
    B=color(3).*ins;
    colorPrior(:,:,1)=colorPrior(:,:,1)+R;
    colorPrior(:,:,2)=colorPrior(:,:,2)+G;
    colorPrior(:,:,3)=colorPrior(:,:,3)+B;
end
R=colorPrior(:,:,1); R=fliplr(R/max(R(:)));
G=colorPrior(:,:,2); G=fliplr(G/max(G(:)));
B=colorPrior(:,:,3); B=fliplr(B/max(B(:)));
colorPrior(:,:,1)=R;
colorPrior(:,:,2)=G;
colorPrior(:,:,3)=B;



  [m,n]=size(R);
%     red = [128 128 128 0 128 192 0 199 0];
%     green = [128 128 64 128 0 128 128 21 0];
%     blue = [128 0 128 0 0 0 128 133 128];
    im=double(sourceImage);
    rgb = zeros(m,n,3);
%     for i = 1:m
%         for j = 1:n
%             if (labelData(i,j) >= 0)
%                 rgb(i,j,1) = 0.5*(red(labelData(i,j)+1))+0.5*im(i,j,1)/255.0;
%                 rgb(i,j,2) = 0.5*(green(labelData(i,j)+1)+im(i,j,2))/255.0;
%                 rgb(i,j,3) = 0.5*(blue(labelData(i,j)+1)+im(i,j,3))/255.0;
%             end
%         end
%     end
% %     figure;
% %     imagesc(rgb);

rgb(:,:,1)=0.5*R+0.5*im(:,:,1)/255.0;
rgb(:,:,2)=0.5*G+0.5*im(:,:,2)/255.0;
rgb(:,:,3)=0.5*B+0.5*im(:,:,3)/255.0;




end