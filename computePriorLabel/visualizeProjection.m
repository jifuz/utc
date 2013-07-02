function visualizeProjection(projection,osm,Cam)

ks=keys(projection);

%figure(1)
hold on
%draw bird view
for i=1:projection.Count
    geometry=projection(ks{i});
    
    for ii=1:numel(geometry.birdView)
        xy=geometry.birdView{ii};
        x=xy(:,1);
        y=xy(:,2);
        plot(x,y);
    end
end

%draw camera
G=Cam.G/(norm(Cam.G));
G=G*100;
% Gl=rotateVectorClockwise(G(1:2),Cam.fov/2);
% Gr=rotateVectorClockwise(G(1:2),-Cam.fov/2);
Gl=(Cam.Gl/norm(Cam.Gl))*100;
Gr=(Cam.Gr/norm(Cam.Gr))*100;
plot([0 G(1)],[0 G(2)],'--k');
plot([0 Gl(1)],[0 Gl(2)],'k');
plot([0 Gr(1)],[0 Gr(2)],'k');


figure
hold on
axis([0 Cam.nx 0 Cam.ny])
%draw projected view
for i=1:projection.Count
    geometry=projection(ks{i});
        
    if isKey(osm.ways,ks{i})
        type=classifyType(osm.ways(ks{i}).tag);
    elseif isKey(osm.nodes,ks{i})
        type=classifyType(osm.nodes(ks{i}).tag);
    else
        error('wrong ID %s',ks{i});
    end
    
    xpyp=geometry.cameraView{1};
    xp=xpyp(:,1); yp=xpyp(:,2);
    
%     osm.ways(ks{i}).tag;
    
    if strcmp(type,'building')
        %get geometry camera view
        xpyptop=geometry.cameraView{2};
        xptop=xpyptop(:,1);
        yptop=xpyptop(:,2);
        
        %filling building color
        for ii=1:min(length(xp),length(xptop))-1 %length(xp)-1 before!!
            polyX=[xptop(ii) xptop(ii+1) xp(ii+1) xp(ii)];
            polyY=[yptop(ii) yptop(ii+1) yp(ii+1) yp(ii)];
            fill(polyX,polyY,[1 0.5 0]);
        end
        
        %plot geometry camera view
        plot(xp,yp,'r','LineWidth',2);
        plot(xptop,yptop,'r','LineWidth',2);
        for ii=1:min(length(xp),length(xptop)) %length(xp) before!!
            plot([xp(ii) xptop(ii)],[yp(ii) yptop(ii)],'r','LineWidth',2);
        end
        
        
    elseif strcmp(type,'highway')
        %get geometry camera view
        xpypl=geometry.cameraView{2};
        xpypr=geometry.cameraView{3};
        xpl=xpypl(:,1);
        ypl=xpypl(:,2);
        xpr=xpypr(:,1);
        ypr=xpypr(:,2);
        
        
        %filling highway color
        for ii=1:min(length(xpl),length(xpr))-1
            polyX=[xpl(ii) xpl(ii+1) xpr(ii+1) xpr(ii) xpl(ii)];
            polyY=[ypl(ii) ypl(ii+1) ypr(ii+1) ypr(ii) ypl(ii)];
            fill(polyX,polyY,[0.5 0.5 0.5]);
        end
        
        %plot geometry camera view
        plot(xp,yp,'k--','LineWidth',1);
        plot(xpl,ypl,'k','LineWidth',1);
        plot(xpr,ypr,'k','LineWidth',1);
        
        %text(mean(xp),mean(yp),[ks{i}]);
        
    elseif strcmp(type,'parking')
        plot(xp,yp,'y','LineWidth',2);
        
    elseif strcmp(type,'bus_stop')
        plot(xp,yp,'b','LineWidth',2);    
    elseif strcmp(type,'traffic_signals')
        plot(xp,yp,'b','LineWidth',2);   
        
    elseif strcmp(type,'undefined')
        %        drawPolyline(x,y);
        %fprintf('shoule never reach here - type:%s\n',type);
        %plot(xp,yp);
    elseif strcmp(type,'park')
        plot(xp,yp,'g','LineWidth',2); 
    else
        %        drawPolyline(x,y);
        fprintf('shoule never reach here - type:%s\n',type);
        %plot(xp,yp);
    end
    
    
end