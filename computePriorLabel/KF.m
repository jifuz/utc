function [KF_xy,KF_dxdy, KF_V]=KF(measurements, model)

global KF_Qnoise
global KF_Rnoise

%parse input
x=measurements.x;
y=measurements.y;
yaw=measurements.yaw;
dx=measurements.dx;
dy=measurements.dy;
dyaw=measurements.dyaw;
t=measurements.time;
t_max=length(x);
start=2;

%model
% Fmodel=[1 0 0 1 0 0;...
%         0 1 0 0 1 0;...
%         0 0 1 0 0 1;...
%         0 0 0 1 0 0;...
%         0 0 0 0 1 0;...
%         0 0 0 0 0 1];
Fmodel=model;

Gmodel=[1 0 0 0;...
        0 1 0 0;];
    
Q=eye(4).*KF_Qnoise;
R=eye(2).*KF_Rnoise;
% R=[1 0;...
%    0 1];


%initial M and V
% M=[x(start); y(start); yaw(start); dx(start); dy(start); dyaw(start)];
M=[x(start); y(start);dx(start); dy(start)];
Mprev=[x(start-1); y(start-1);dx(start-1); dy(start-1)];

V=eye(4);
MM=cell(t_max,1);
VV=cell(t_max,1);
MM{1}=M;
VV{1}=V;

%KF main loop
for i=start:t_max
    
   %predict
   Mprev=M;
    [M,V]=ekf_predict(M,Mprev,V,Fmodel,Q,i);
    
    %update
    [M,V]=ekf_update(M,V,Gmodel,measurements,R,i);
    
    MM{i}=M';
    VV{i}=V;
    
end


KF_xy=MM(start:end,1);
KF_xy=cell2mat(KF_xy);
KF_dxdy=KF_xy(:,3:4);
KF_xy=KF_xy(:,1:2);
KF_V=VV;
end



function [M_,V_]=ekf_predict(M,Mprev,V,Fmodel,Q,i)

%M
%M_=Fmodel*M;
M_=driveOnNearestRoadSeg(Fmodel,M,Mprev,i);

%V
Jf=computeJf(M_,M);
dF=Jf;
V_=dF*V*dF'+Q;

end

 
function Jf=computeJf(f,x)

% m=length(f);
% n=length(x);
% Jf=zeros(m,n);
% for i=1:m
%     for j=1:n
%         Jf=(f(i)-x(i))
%     end
% end

Jf=f/x;

end


function [M_,V_]=ekf_update(M,V,Gmodel,measurements,R,i)

% Gm=driveOnNearestRoadSeg(Gmodel,M);
C=Gmodel;
Gm=(C*M);
mx=measurements.x;
my=measurements.y;
measurement=[mx(i);my(i)];
innovation=(measurement-Gm);

Jg=C;


gain=V*Jg' / (Jg*V*Jg'+R);
M_=M+gain*innovation;
V_=(eye(4)-gain*Jg)*V;

end

function M_=driveOnNearestRoadSeg(Fmodel,M,Mprev,i)

x=M(1);
y=M(2);
dx=M(3);
dy=M(4);
% find the nearest lane seg
minD=9999;
lanes=Fmodel.laneCell;
p=[x y];
for i=1:length(lanes)
    lane=lanes{i};
    for j=1:length(lane)-1
        a=lane(j,:);
        b=lane(j+1,:);
        d=distancePointToLineSegment(p, a, b);
        if d<minD
            minInd=[i,j];
            mina=a;
            minb=b;
        end
    end
end


yaw=[dx dy];

v1=(mina-minb);
v2=-v1;
D=pdist2(yaw, [v1;v2]);
if D(1)>D(2)
    v=v2;
else
    v=v1;
end

%project yaw on v
yaw_onV=(dot(yaw,v)/norm(v)^2)*v;


%compute predicted x y dx dy
deltaT=Fmodel.time(i);

x_=x+deltaT*yaw_onV(1);
y_=y+deltaT*yaw_onV(2);

% dx_=mean([(x-Mprev(1))/deltaT, yaw_onV(1)]);
% dy_=mean([(y-Mprev(2))/deltaT, yaw_onV(2)]);
dx_=yaw_onV(1);
dy_=yaw_onV(2);

M_=[x_; y_; dx_; dy_];
end



function distance = distancePointToLineSegment(p, a, b)
    % p is a point
    % a is the start of a line segment
    % b is the end of a line segment

    pa_distance = sum((a-p).^2);
    pb_distance = sum((b-p).^2);

    unitab = (a-b) / norm(a-b);

    if pa_distance < pb_distance
        d = dot((p-a)/norm(p-a), -unitab);

        if d > 0
            distance = sqrt(pa_distance * (1 - d*d));
        else
            distance = sqrt(pa_distance);
        end
    else
        d = dot((p-b)/norm(p-b), unitab);

        if d > 0
            distance = sqrt(pb_distance * (1 - d*d));
        else
            distance = sqrt(pb_distance);
        end
    end
end