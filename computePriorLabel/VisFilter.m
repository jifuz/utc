function [xVis,yVis,ind]=VisFilter(Cam,x,y)
%
% figure
% hold on
% plot(x,y);
% plot(0,0,'x');

xVis=x;
yVis=y;

nx=Cam.nx;
ny=Cam.ny;
n=Cam.n;
fov=Cam.fov;
f=Cam.f;
E=Cam.E;
G=Cam.G;
T=Cam.T;

ind=x.*0;
vertGL=vertVectorL(G);
vertGR=vertVectorR(G);

Gl=Cam.Gl;
Gr=Cam.Gr;
for i=1:length(x)
    [xg,yg]=projectOnVector([x(i),y(i)],G(1:2));
    if G(1)/xg<0 || G(2)/yg<0 || (xg==0 && yg==0) %if behind in screen
        [xvg,yvg]=projectOnVector([x(i),y(i)],vertGR);%project on screen along x-axis of camera
        
        %%%%%%%
%         [xgl,ygl]=projectOnVector([x(i),y(i)],vertGL);
%         [xgr,ygr]=projectOnVector([x(i),y(i)],vertGR);
%         
%         %leftGaze-projected is negative, use project on right instead
%         if Gl(1)/xgl<0 || Gl(2)/ygl<0 || (xgl==0 && ygl==0)
%             [xgr,ygr]=projectOnVector([x(i),y(i)],Gr);
%             xVis(i)=xgr;
%             yVis(i)=ygr;
%             
%         %rightGaze-projected is negative, use project on left G instead
%         elseif Gr(1)/xgr<0 || Gr(2)/ygr<0 || (xgr==0 && ygr==0)
%             [xgl,ygl]=projectOnVector([x(i),y(i)],Gl);
%             xVis(i)=xgl;
%             yVis(i)=ygl;
%         end
        %%%%%
        
        %adjust to project on camera screen a bit positive on y-axis
        %adj=G(1:2)-[xvg,yvg];
        adj=G(1:2);
        adj=adj/norm(adj);
        xvg=xvg+0.01*adj(1);
        yvg=yvg+0.01*adj(2);
        
        xVis(i)=xvg;
        yVis(i)=yvg;
        ind(i)=1;
    end
    
end

% figure
% hold on
% plot(xVis,yVis);
% plot(0,0,'x');

end


%project A onto B
function [C1,C2]=projectOnVector(A,B)
C = (dot(A,B)/norm(B)^2)*B;
C1=C(1);
C2=C(2);
end

%[0 1] -> [-1 0]
%[1 1] -> [-1 1]
function vertG=vertVectorL(G)
vertG=-[G(2) -G(1)];
end

function vertG=vertVectorR(G)
vertG=-vertVectorL(G);
end