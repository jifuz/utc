function [labelLayers,layerPrior]=visualizeProjection_fillingLabels(reprojection,osm,Cam)

global numSample

%% area layer, road layer, building layer
labelLayers=containers.Map;
ks=keys(reprojection);

%layers={'parking','park','highway','footway','building'};
layers={'building','highway'};

keysLayers=containers.Map;
keysLayers('parking')={};
keysLayers('park')={};
keysLayers('highway')={};
keysLayers('footway')={};
keysLayers('building')={};


for i=1:reprojection.Count
    %     geometry=reprojection(ks{i});
    
    if ~isKey(osm.ways,ks{i})
        continue;
    end
    type=classifyType(osm.ways(ks{i}).tag);
    
    currCell=keysLayers(type);
    currCell{end+1}=ks{i};
    keysLayers(type)=currCell;
    
    %     if strcmp(type,'parking')
    %         parkingKeys{end+1}=ks{i};
    %     elseif strcmp(type,'park')
    %         parkKeys{end+1}=ks{i};
    %
    %     elseif strcmp(type,'highway')
    %         highwayKeys{end+1}=ks{i};
    %
    %     elseif strcmp(type,'footway')
    %         footwayKeys{end+1}=ks{i};
    %
    %     elseif strcmp(type,'building')
    %         buildingKeys{end+1}=ks{i};
    %
    %     else
    %
    %     end
    
    %     keys_sort={parkingKeys,parkKeys,highwayKeys,footwayKeys,buildingKeys};
    
end


for i=1:length(layers)
    currLayer=layers{i};
    labelLayers(currLayer)=zeros(240,795);
end

%% highway layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% ins=cell(1,numSample);
ins=zeros(240,795);
% Y=meshgrid(1:Cam.nx,1:Cam.ny);
% X=meshgrid(1:Cam.ny,1:Cam.nx);X=X';
Y=meshgrid(1:Cam.ny,1:Cam.nx);
X=meshgrid(1:Cam.nx,1:Cam.ny);X=X';


for l=1:length(layers)
    fprintf('%s, %d/%d...\n',currLayer,l,length(layers));
    ins=zeros(240,795);
    currLayer=layers{l};
    currLayerKeys=keysLayers(currLayer);
    parfor iii=1:numSample
        fprintf('KF:%d/%d...',iii,numSample);
        currIn=logical(zeros(240,795));
        
        for i=1:length(currLayerKeys)
            %                 fprintf('%s, %d/%d...\n',currLayer,i,length(currLayerKeys));
            geometry=reprojection(currLayerKeys{i});
            xpyps=geometry.cameraView{1};
            xp=xpyps{iii,1};
            yp=xpyps{iii,2};
            
            
            if strcmp(currLayer,'highway') || strcmp(currLayer,'footway')
                %get geometry camera view
                xpypls=geometry.cameraView{2};
                xpyprs=geometry.cameraView{3};
                xpl=xpypls{iii,1};
                ypl=xpypls{iii,2};
                xpr=xpyprs{iii,1};
                ypr=xpyprs{iii,2};
                
                for ii=1:length(xpl)-1
                    polyX=[xpl(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
                    polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
                    in=inpolygon(X,Y,[polyX polyX(1)],[polyY polyY(1)]);
                    in_=in';
                    in_=imresize(in_,[240 795]);
                    currIn=currIn | in_;
                end
                
                
            elseif strcmp(currLayer,'building')
                xpypstop=geometry.cameraView{2};
                xptop=xpypstop{iii,1};
                yptop=xpypstop{iii,2};
                
                for ii=1:length(xp)-1
                    polyX=[xptop(ii) xptop(ii+1) xp(ii+1) xp(ii)];
                    polyY=[yptop(ii) yptop(ii+1) yp(ii+1) yp(ii)];
                    in=inpolygon(X,Y,[polyX polyX(1)],[polyY polyY(1)]);
                    %                     in=inpolygon(Y,X,[polyY polyY(1)],[polyX polyX(1)]);
                    in_=in';
                    in_=imresize(in_,[240 795]);
                    currIn=currIn | in_;
                end
                
                
            else
                polyX=xp;
                polyY=yp;
                in=inpolygon(X,Y,[polyX polyX(1)],[polyY polyY(1)]);
                in_=in';
                in_=imresize(in_,[240 795]);
                currIn=currIn | in_;
            end
            
            
        end
        currIn=flipud(currIn);
        ins=ins+double(currIn);
        
        
        %         ins_=flipud(ins');ins_=ins_./max(ins_(:));ins_(ins_<0.5)=0;ins_=exp(ins_);imagesc(ins_);
        %         figure
        %         axis([0 Cam.nx 0 Cam.ny])
        %         subplot(121)
        %         imagesc(flipud(ins'))
        %         subplot(122)
        %         imagesc(ins_);
        %         title(currLayer)
%         figure
%         imagesc(currIn)
        
    end
    ins(ins>numSample)=numSample;
    ins_=flipud(ins');ins_=ins_./numSample;
    %ins_=exp(ins_.*10)/exp(10);%%%!!!Test exp distribution
    labelLayers(currLayer)=ins_;
    %     ins=ins+double(currIn);
end


% for l=1:length(layers)
%     currLayer=layers{l};
%     ins=labelLayers(currLayer);
ins_=flipud(ins');ins_=ins_./numSample;ins_(ins_<0.3)=0;ins_=exp(ins_);imagesc(ins_);
%     ins_=flipud(ins');ins_=ins_./max(ins_(:));ins_(ins_<0.3)=0;imagesc(ins_);
%
%     figure
%     axis([0 Cam.ny 0 Cam.nx])
%     %         subplot(121)
%     imagesc(flipud(ins'))
%     %         subplot(122)
%     title(currLayer)
%
%     figure
%     imagesc(ins_);
%     title(currLayer)
% end

buildingLayer=labelLayers('building');
highwayLayer=labelLayers('highway');
highwayLayer(buildingLayer>0.8)=0;
labelLayers('highway')=highwayLayer;

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
    %ins(ins<0.2)=0;
    layerPrior(currLayer)=ins;
    color=layerColor(currLayer);
    R=color(1).*ins;
    G=color(2).*ins;
    B=color(3).*ins;
    colorPrior(:,:,1)=colorPrior(:,:,1)+R;
    colorPrior(:,:,2)=colorPrior(:,:,2)+G;
    colorPrior(:,:,3)=colorPrior(:,:,3)+B;
end
R=colorPrior(:,:,1);
G=colorPrior(:,:,2);
B=colorPrior(:,:,3);
colorPrior(:,:,1)=fliplr(R/max(R(:)));
colorPrior(:,:,2)=fliplr(G/max(G(:)));
colorPrior(:,:,3)=fliplr(B/max(B(:)));

figure
imagesc(colorPrior)
axis equal

return
% 
% for l=1:length(layers)
%     currLayer=layers{l};
%     ins=labelLayers(currLayer);
%     ins_=flipud(ins');ins_=ins_./max(ins_(:));ins_(ins_<0.5)=0;ins_=exp(ins_);imagesc(ins_);
%     figure
%     axis([0 Cam.nx 0 Cam.ny])
%     subplot(121)
%     imagesc(flipud(ins'))
%     subplot(122)
%     imagesc(ins_);
%     title(currLayer)
% end
% 
% 
% 
% %convert to label layer
% ins_=flipud(ins');ins_=ins_./max(ins_(:));ins_(ins_<0.5)=0;ins_=exp(ins_);imagesc(ins_);
% figure
% imagesc(flipud(ins'))
% return
% 
% 
% 
% 
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%     
%     if ~isKey(osm.ways,ks{i})
%         continue;
%     end
%     type=classifyType(osm.ways(ks{i}).tag);
%     
%     xpyp=geometry.cameraView{1};
%     xp=xpyp(:,1); yp=xpyp(:,2);
%     
%     osm.ways(ks{i}).tag
%     
%     if strcmp(type,'building')
%         %get geometry camera view
%         xpyptop=geometry.cameraView{2};
%         xptop=xpyptop(:,1);
%         yptop=xpyptop(:,2);
%         
%         
%         %filling color
%         for ii=1:length(xp)-1
%             polyX=[xptop(ii) xptop(ii+1) xp(ii+1) xp(ii)];
%             polyY=[yptop(ii) yptop(ii+1) yp(ii+1) yp(ii)];
%             fill(polyX,polyY,[0 0 0]);
%             filled=1;
%         end
%         
%         
%         %plot geometry camera view
%         %         plot(xp,yp,'r','LineWidth',2);
%         %         plot(xptop,yptop,'r','LineWidth',2);
%         %         for ii=1:length(xp)
%         %             plot([xp(ii) xptop(ii)],[yp(ii) yptop(ii)],'r','LineWidth',2);
%         %         end
%         
%         
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %draw projected view and fill labels
% %% park layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% filled=0;
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%     
%     if ~isKey(osm.ways,ks{i})
%         continue;
%     end
%     type=classifyType(osm.ways(ks{i}).tag);
%     
%     xpyp=geometry.cameraView{1};
%     xp=xpyp(:,1); yp=xpyp(:,2);
%     
%     osm.ways(ks{i}).tag
%     
%     if strcmp(type,'park')
%         
%         %filling color
%         fill(xp,yp,[0 1 0]);
%         filled=1;
%         in=inpolygon(Cam.nx,Cam.ny,xp,yp);
%         %
%         %         for ii=1:length(xp)-1
%         %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
%         %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
%         %         end
%         
%     end
%     
% end
% 
% %convert to label layer
% if filled
%     axis off;
%     F=getframe;
%     im=imresize(F.cdata, [Cam.nx Cam.ny]);
%     imgray=1-(double(rgb2gray(im))./255);
%     imgray(imgray<0.89)=0;
%     labelLayers('park')=imgray;
%     %     figure
%     %     imshow(imgray)
% end
% close
% 
% 
% %% area layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% filled=0
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%     
%     if ~isKey(osm.ways,ks{i})
%         continue;
%     end
%     type=classifyType(osm.ways(ks{i}).tag);
%     
%     xpyp=geometry.cameraView{1};
%     xp=xpyp(:,1); yp=xpyp(:,2);
%     
%     osm.ways(ks{i}).tag
%     
%     if strcmp(type,'parking')
%         
%         %filling color
%         fill(xp,yp,[0 0 0]);
%         filled=1;
%         in=inpolygon(Cam.nx,Cam.ny,xp,yp);
%         
%         %
%         %         for ii=1:length(xp)-1
%         %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
%         %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
%         %         end
%         
%     end
%     
% end
% 
% %convert to label layer
% if filled
%     axis off;
%     F=getframe;
%     im=imresize(F.cdata, [Cam.nx Cam.ny]);
%     imgray=1-(double(rgb2gray(im))./255);
%     imgray(imgray<0.89)=0;
%     labelLayers('parking')=imgray;
%     %     figure
%     %     imshow(imgray)
% end
% close
% 
% %% highway layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% filled=0;
% ins=cell(1,numSample);
% ins=zeros(240,795);
% for iii=1:numSample
%     currIn=logical(zeros(240,795));
%     fprintf('%d/%d...\n',iii,numSample);
%     for i=1:reprojection.Count
%         %     fprintf('%d/%d...\n',i,reprojection.Count);
%         geometry=reprojection(ks{i});
%         
%         if ~isKey(osm.ways,ks{i})
%             continue;
%         end
%         type=classifyType(osm.ways(ks{i}).tag);
%         
%         xpyps=geometry.cameraView{1};
%         %iii=2;
%         %     inCurr=[];
%         %     for iii=1:length(xpyps)
%         xp=xpyps{iii,1};
%         yp=xpyps{iii,2};
%         
%         %          xpyp=geometry.cameraView{1};
%         %          xp=xpyp(:,1); yp=xpyp(:,2);
%         
%         %         osm.ways(ks{i}).tag
%         
%         if strcmp(type,'highway')
%             %get geometry camera view
%             xpypls=geometry.cameraView{2};
%             xpyprs=geometry.cameraView{3};
%             %xpypl=xpypls{iii};
%             xpl=xpypls{iii,1};
%             ypl=xpypls{iii,2};
%             xpr=xpyprs{iii,1};
%             ypr=xpyprs{iii,2};
%             
%             %             xpypl=geometry.cameraView{2};
%             %         xpypr=geometry.cameraView{3};
%             %         xpl=xpypl(:,1);
%             %         ypl=xpypl(:,2);
%             %         xpr=xpypr(:,1);
%             %         ypr=xpypr(:,2);
%             %
%             
%             %filling color
%             
%             for ii=1:length(xpl)-1
%                 polyX=[xpl(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
%                 polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
%                 fill(polyX,polyY,[0 0 0]);
%                 filled=1;
%                 %             polyX(polyX>Cam.nx)=Cam.nx;
%                 %             polyY(polyY>Cam.ny)=Cam.ny;
%                 %             polyX(polyX<0)=0;
%                 %             polyY(polyY<0)=0;
%                 
%                 Y=meshgrid(1:Cam.ny,1:Cam.nx);
%                 X=meshgrid(1:Cam.nx,1:Cam.ny);X=X';
%                 
%                 in=inpolygon(X,Y,[polyX polyX(1)],[polyY polyY(1)]);
%                 
%                 %                 if length(inCurr)<1;
%                 %                 inCurr=double(in);
%                 %                 else
%                 %                     inCurr=double(in)+inCurr;
%                 %                 end
%                 
%                 %             polyX=polyX';
%                 %             polyY=polyY';
%                 %             X=X';
%                 %             Y=Y';
%                 %             in=in';
%                 
%                 %             plot(polyX,polyY,X(in),Y(in),'r+',X(~in),Y(~in),'b.')
%                 %                 plot(polyX,polyY,X(in),Y(in),'r+')
%                 
%                 currIn=currIn | in;
%                 
%                 %                 max(in(:))
%             end
%             
%             %plot geometry camera view
%             %          plot(xp,yp,'k--');
%             %          plot(xpl,ypl,'k');
%             %          plot(xpr,ypr,'k');
%             %
%             
%             
%             %text(mean(xp),mean(yp),[ks{i}]);
%             
%             %     elseif strcmp(type,'parking')
%             %         plot(xp,yp,'y','LineWidth',2);
%             %
%             %     elseif strcmp(type,'undefined')
%             %         %        drawPolyline(x,y);
%             %         %fprintf('shoule never reach here - type:%s\n',type);
%             %         %plot(xp,yp);
%         else
%             %        drawPolyline(x,y);
%             fprintf('shoule never reach here - type:%s\n',type);
%             %plot(xp,yp);
%         end
%         
%         %     end
%         %     ins{i}=inCurr;
%     end
%     % imshow(double(currIn))
%     % hold off
%     % plot(0,0,'b')
%     % hold on
%     ins=ins+double(currIn);
% end
% %convert to label layer
% ins_=flipud(ins');ins_=ins_./max(ins_(:));ins_(ins_<0.5)=0;ins_=exp(ins_);imagesc(ins_);
% figure
% imagesc(flipud(ins'))
% return
% if filled
%     axis off;
%     F=getframe;
%     im=imresize(F.cdata, [Cam.nx Cam.ny]);
%     imgray=1-(double(rgb2gray(im))./255);
%     imgray(imgray<0.89)=0;
%     labelLayers('highway')=imgray;
%     %     figure
%     %     imshow(imgray)
% end
% close
% 
% %% footway layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% filled=0;
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%     
%     if ~isKey(osm.ways,ks{i})
%         continue;
%     end
%     type=classifyType(osm.ways(ks{i}).tag);
%     
%     xpyp=geometry.cameraView{1};
%     xp=xpyp(:,1); yp=xpyp(:,2);
%     
%     osm.ways(ks{i}).tag
%     
%     if strcmp(type,'footway')
%         %get geometry camera view
%         xpypl=geometry.cameraView{2};
%         xpypr=geometry.cameraView{3};
%         xpl=xpypl(:,1);
%         ypl=xpypl(:,2);
%         xpr=xpypr(:,1);
%         ypr=xpypr(:,2);
%         
%         %filling color
%         for ii=1:length(xpl)-1
%             polyX=[xpl(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
%             polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
%             fill(polyX,polyY,[0 0 0]);
%             filled=1;
%         end
%         
%         %plot geometry camera view
%         %          plot(xp,yp,'k--');
%         %          plot(xpl,ypl,'k');
%         %          plot(xpr,ypr,'k');
%         %
%         
%         
%         %text(mean(xp),mean(yp),[ks{i}]);
%         
%         %     elseif strcmp(type,'parking')
%         %         plot(xp,yp,'y','LineWidth',2);
%         %
%         %     elseif strcmp(type,'undefined')
%         %         %        drawPolyline(x,y);
%         %         %fprintf('shoule never reach here - type:%s\n',type);
%         %         %plot(xp,yp);
%     else
%         %        drawPolyline(x,y);
%         fprintf('shoule never reach here - type:%s\n',type);
%         %plot(xp,yp);
%     end
%     
% end
% 
% %convert to label layer
% if filled
%     axis off;
%     F=getframe;
%     im=imresize(F.cdata, [Cam.nx Cam.ny]);
%     imgray=1-(double(rgb2gray(im))./255);
%     imgray(imgray<0.89)=0;
%     labelLayers('footway')=imgray;
%     %     figure
%     %     imshow(imgray)
% end
% close
% 
% %% building layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% filled=0;
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%     
%     if ~isKey(osm.ways,ks{i})
%         continue;
%     end
%     type=classifyType(osm.ways(ks{i}).tag);
%     
%     xpyp=geometry.cameraView{1};
%     xp=xpyp(:,1); yp=xpyp(:,2);
%     
%     osm.ways(ks{i}).tag
%     
%     if strcmp(type,'building')
%         %get geometry camera view
%         xpyptop=geometry.cameraView{2};
%         xptop=xpyptop(:,1);
%         yptop=xpyptop(:,2);
%         
%         
%         %filling color
%         for ii=1:length(xp)-1
%             polyX=[xptop(ii) xptop(ii+1) xp(ii+1) xp(ii)];
%             polyY=[yptop(ii) yptop(ii+1) yp(ii+1) yp(ii)];
%             fill(polyX,polyY,[0 0 0]);
%             filled=1;
%         end
%         
%         
%         %plot geometry camera view
%         %         plot(xp,yp,'r','LineWidth',2);
%         %         plot(xptop,yptop,'r','LineWidth',2);
%         %         for ii=1:length(xp)
%         %             plot([xp(ii) xptop(ii)],[yp(ii) yptop(ii)],'r','LineWidth',2);
%         %         end
%         
%         
%     end
% end
% 
% %convert to label layer
% if filled
%     axis off;
%     F=getframe;
%     im=imresize(F.cdata, [Cam.nx Cam.ny]);
%     imgray=1-(double(rgb2gray(im))./255);
%     imgray(imgray<0.89)=0;
%     labelLayers('building')=imgray;
%     %     figure
%     %     imshow(imgray)
% end
% close
% 
% 
% 
% %% node layer
% figure
% hold on
% axis([0 Cam.nx 0 Cam.ny])
% filled=0;
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%     
%     
%     if ~isKey(osm.nodes,ks{i})
%         continue;
%     end
%     type=classifyType(osm.nodes(ks{i}).tag);
%     
%     xpyp=geometry.cameraView{1};
%     xp=xpyp(:,1); yp=xpyp(:,2);
%     
%     osm.nodes(ks{i}).tag
%     
%     if strcmp(type,'bus_stop')
%         
%         %filling color
%         fill(xp,yp,[0.5 0 0.5]);
%         filled=1;
%         %
%         %         for ii=1:length(xp)-1
%         %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
%         %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
%         %         end
%     elseif strcmp(type,'traffic_signals')
%         
%         %filling color
%         fill(xp,yp,[0.5 0 0.5]);
%         filled=1;
%         %
%         %         for ii=1:length(xp)-1
%         %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
%         %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
%         %
%     end
% end
% 
% %convert to label layer
% if filled
%     axis off;
%     F=getframe;
%     im=imresize(F.cdata, [Cam.nx Cam.ny]);
%     imgray=1-((rgb2gray(im))./255);
%     imgray(imgray<0.89)=0;
%     labelLayers('transportation_objects')=imgray;
%     %          figure
%     %          imshow(imgray)
% end
% close


end


