function [labelLayers]=visualizeProjection_fillingLabels(reprojection,osm,Cam)

numSample=100;

%% area layer, road layer, building layer
labelLayers=containers.Map;
ks=keys(reprojection);

% figure
% hold on
% axis equal
% %draw bird view
% for i=1:reprojection.Count
%     geometry=reprojection(ks{i});
%
%     for ii=1:numel(geometry.birdView)
%         xy=geometry.birdView{ii};
%         x=xy(:,1);
%         y=xy(:,2);
%         %[x,y]=rotateToNorth(x,y,Cam);
%         plot(x,y);
%     end
% end
% %draw camera
% G=Cam.G/(norm(Cam.G));
% G=G*100;
% Gl=rotateVectorClockwise(G(1:2),Cam.fov/2);
% Gr=rotateVectorClockwise(G(1:2),-Cam.fov/2);
% plot([0 G(1)],[0 G(2)],'--g');
% plot([0 Gl(1)],[0 Gl(2)],'g');
% plot([0 Gr(1)],[0 Gr(2)],'g');
%





%draw projected view and fill labels
%% park layer
figure
hold on
axis([0 Cam.nx 0 Cam.ny])
filled=0;
for i=1:reprojection.Count
    geometry=reprojection(ks{i});
    
    if ~isKey(osm.ways,ks{i})
        continue;
    end
    type=classifyType(osm.ways(ks{i}).tag);
    
    xpyp=geometry.cameraView{1};
    xp=xpyp(:,1); yp=xpyp(:,2);
    
    osm.ways(ks{i}).tag
    
    if strcmp(type,'park')
        
        %filling color
        fill(xp,yp,[0 1 0]);
        filled=1;
        in=inpolygon(Cam.nx,Cam.ny,xp,yp);
        %
        %         for ii=1:length(xp)-1
        %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
        %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
        %         end
        
    end
    
end

%convert to label layer
if filled
    axis off;
    F=getframe;
    im=imresize(F.cdata, [Cam.nx Cam.ny]);
    imgray=1-(double(rgb2gray(im))./255);
    imgray(imgray<0.89)=0;
    labelLayers('park')=imgray;
    %     figure
    %     imshow(imgray)
end
close


%% area layer
figure
hold on
axis([0 Cam.nx 0 Cam.ny])
filled=0
for i=1:reprojection.Count
    geometry=reprojection(ks{i});
    
    if ~isKey(osm.ways,ks{i})
        continue;
    end
    type=classifyType(osm.ways(ks{i}).tag);
    
    xpyp=geometry.cameraView{1};
    xp=xpyp(:,1); yp=xpyp(:,2);
    
    osm.ways(ks{i}).tag
    
    if strcmp(type,'parking')
        
        %filling color
        fill(xp,yp,[0 0 0]);
        filled=1;
        in=inpolygon(Cam.nx,Cam.ny,xp,yp);
        
        %
        %         for ii=1:length(xp)-1
        %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
        %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
        %         end
        
    end
    
end

%convert to label layer
if filled
    axis off;
    F=getframe;
    im=imresize(F.cdata, [Cam.nx Cam.ny]);
    imgray=1-(double(rgb2gray(im))./255);
    imgray(imgray<0.89)=0;
    labelLayers('parking')=imgray;
    %     figure
    %     imshow(imgray)
end
close

%% highway layer
figure
hold on
axis([0 Cam.nx 0 Cam.ny])
filled=0;
ins=cell(1,numSample);
ins=zeros(240,795);
for iii=1:numSample
    currIn=logical(zeros(240,795));
    fprintf('%d/%d...\n',iii,numSample);
    for i=1:reprojection.Count
        %     fprintf('%d/%d...\n',i,reprojection.Count);
        geometry=reprojection(ks{i});
        
        if ~isKey(osm.ways,ks{i})
            continue;
        end
        type=classifyType(osm.ways(ks{i}).tag);
        
        xpyps=geometry.cameraView{1};
        %iii=2;
        %     inCurr=[];
        %     for iii=1:length(xpyps)
        xp=xpyps{iii,1};
        yp=xpyps{iii,2};
        
        %          xpyp=geometry.cameraView{1};
        %          xp=xpyp(:,1); yp=xpyp(:,2);
        
        %         osm.ways(ks{i}).tag
        
        if strcmp(type,'highway')
            %get geometry camera view
            xpypls=geometry.cameraView{2};
            xpyprs=geometry.cameraView{3};
            %xpypl=xpypls{iii};
            xpl=xpypls{iii,1};
            ypl=xpypls{iii,2};
            xpr=xpyprs{iii,1};
            ypr=xpyprs{iii,2};
            
            %             xpypl=geometry.cameraView{2};
            %         xpypr=geometry.cameraView{3};
            %         xpl=xpypl(:,1);
            %         ypl=xpypl(:,2);
            %         xpr=xpypr(:,1);
            %         ypr=xpypr(:,2);
            %
            
            %filling color
            
            for ii=1:length(xpl)-1
                polyX=[xpl(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
                polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
                fill(polyX,polyY,[0 0 0]);
                filled=1;
                %             polyX(polyX>Cam.nx)=Cam.nx;
                %             polyY(polyY>Cam.ny)=Cam.ny;
                %             polyX(polyX<0)=0;
                %             polyY(polyY<0)=0;
                
                Y=meshgrid(1:Cam.ny,1:Cam.nx);
                X=meshgrid(1:Cam.nx,1:Cam.ny);X=X';
                
                in=inpolygon(X,Y,[polyX polyX(1)],[polyY polyY(1)]);
                
                %                 if length(inCurr)<1;
                %                 inCurr=double(in);
                %                 else
                %                     inCurr=double(in)+inCurr;
                %                 end
                
                %             polyX=polyX';
                %             polyY=polyY';
                %             X=X';
                %             Y=Y';
                %             in=in';
                
                %             plot(polyX,polyY,X(in),Y(in),'r+',X(~in),Y(~in),'b.')
                %                 plot(polyX,polyY,X(in),Y(in),'r+')
                
                currIn=currIn | in;
                
                %                 max(in(:))
            end
            
            %plot geometry camera view
            %          plot(xp,yp,'k--');
            %          plot(xpl,ypl,'k');
            %          plot(xpr,ypr,'k');
            %
            
            
            %text(mean(xp),mean(yp),[ks{i}]);
            
            %     elseif strcmp(type,'parking')
            %         plot(xp,yp,'y','LineWidth',2);
            %
            %     elseif strcmp(type,'undefined')
            %         %        drawPolyline(x,y);
            %         %fprintf('shoule never reach here - type:%s\n',type);
            %         %plot(xp,yp);
        else
            %        drawPolyline(x,y);
            fprintf('shoule never reach here - type:%s\n',type);
            %plot(xp,yp);
        end
        
        %     end
        %     ins{i}=inCurr;
    end
    % imshow(double(currIn))
    % hold off
    % plot(0,0,'b')
    % hold on
    ins=ins+double(currIn);
end
%convert to label layer
ins_=flipud(ins');ins_=ins_./max(ins_(:));ins_(ins_<0.5)=0;ins_=exp(ins_);imagesc(ins_);
figure
imagesc(flipud(ins'))
return
if filled
    axis off;
    F=getframe;
    im=imresize(F.cdata, [Cam.nx Cam.ny]);
    imgray=1-(double(rgb2gray(im))./255);
    imgray(imgray<0.89)=0;
    labelLayers('highway')=imgray;
    %     figure
    %     imshow(imgray)
end
close

%% footway layer
figure
hold on
axis([0 Cam.nx 0 Cam.ny])
filled=0;
for i=1:reprojection.Count
    geometry=reprojection(ks{i});
    
    if ~isKey(osm.ways,ks{i})
        continue;
    end
    type=classifyType(osm.ways(ks{i}).tag);
    
    xpyp=geometry.cameraView{1};
    xp=xpyp(:,1); yp=xpyp(:,2);
    
    osm.ways(ks{i}).tag
    
    if strcmp(type,'footway')
        %get geometry camera view
        xpypl=geometry.cameraView{2};
        xpypr=geometry.cameraView{3};
        xpl=xpypl(:,1);
        ypl=xpypl(:,2);
        xpr=xpypr(:,1);
        ypr=xpypr(:,2);
        
        %filling color
        for ii=1:length(xpl)-1
            polyX=[xpl(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
            fill(polyX,polyY,[0 0 0]);
            filled=1;
        end
        
        %plot geometry camera view
        %          plot(xp,yp,'k--');
        %          plot(xpl,ypl,'k');
        %          plot(xpr,ypr,'k');
        %
        
        
        %text(mean(xp),mean(yp),[ks{i}]);
        
        %     elseif strcmp(type,'parking')
        %         plot(xp,yp,'y','LineWidth',2);
        %
        %     elseif strcmp(type,'undefined')
        %         %        drawPolyline(x,y);
        %         %fprintf('shoule never reach here - type:%s\n',type);
        %         %plot(xp,yp);
    else
        %        drawPolyline(x,y);
        fprintf('shoule never reach here - type:%s\n',type);
        %plot(xp,yp);
    end
    
end

%convert to label layer
if filled
    axis off;
    F=getframe;
    im=imresize(F.cdata, [Cam.nx Cam.ny]);
    imgray=1-(double(rgb2gray(im))./255);
    imgray(imgray<0.89)=0;
    labelLayers('footway')=imgray;
    %     figure
    %     imshow(imgray)
end
close

%% building layer
figure
hold on
axis([0 Cam.nx 0 Cam.ny])
filled=0;
for i=1:reprojection.Count
    geometry=reprojection(ks{i});
    
    if ~isKey(osm.ways,ks{i})
        continue;
    end
    type=classifyType(osm.ways(ks{i}).tag);
    
    xpyp=geometry.cameraView{1};
    xp=xpyp(:,1); yp=xpyp(:,2);
    
    osm.ways(ks{i}).tag
    
    if strcmp(type,'building')
        %get geometry camera view
        xpyptop=geometry.cameraView{2};
        xptop=xpyptop(:,1);
        yptop=xpyptop(:,2);
        
        
        %filling color
        for ii=1:length(xp)-1
            polyX=[xptop(ii) xptop(ii+1) xp(ii+1) xp(ii)];
            polyY=[yptop(ii) yptop(ii+1) yp(ii+1) yp(ii)];
            fill(polyX,polyY,[0 0 0]);
            filled=1;
        end
        
        
        %plot geometry camera view
        %         plot(xp,yp,'r','LineWidth',2);
        %         plot(xptop,yptop,'r','LineWidth',2);
        %         for ii=1:length(xp)
        %             plot([xp(ii) xptop(ii)],[yp(ii) yptop(ii)],'r','LineWidth',2);
        %         end
        
        
    end
end

%convert to label layer
if filled
    axis off;
    F=getframe;
    im=imresize(F.cdata, [Cam.nx Cam.ny]);
    imgray=1-(double(rgb2gray(im))./255);
    imgray(imgray<0.89)=0;
    labelLayers('building')=imgray;
    %     figure
    %     imshow(imgray)
end
close



%% node layer
figure
hold on
axis([0 Cam.nx 0 Cam.ny])
filled=0;
for i=1:reprojection.Count
    geometry=reprojection(ks{i});
    
    
    if ~isKey(osm.nodes,ks{i})
        continue;
    end
    type=classifyType(osm.nodes(ks{i}).tag);
    
    xpyp=geometry.cameraView{1};
    xp=xpyp(:,1); yp=xpyp(:,2);
    
    osm.nodes(ks{i}).tag
    
    if strcmp(type,'bus_stop')
        
        %filling color
        fill(xp,yp,[0.5 0 0.5]);
        filled=1;
        %
        %         for ii=1:length(xp)-1
        %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
        %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
        %         end
    elseif strcmp(type,'traffic_signals')
        
        %filling color
        fill(xp,yp,[0.5 0 0.5]);
        filled=1;
        %
        %         for ii=1:length(xp)-1
        %            polyX=[xp(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
        %            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
        %
    end
end

%convert to label layer
if filled
    axis off;
    F=getframe;
    im=imresize(F.cdata, [Cam.nx Cam.ny]);
    imgray=1-((rgb2gray(im))./255);
    imgray(imgray<0.89)=0;
    labelLayers('transportation_objects')=imgray;
    %          figure
    %          imshow(imgray)
end
close


end


