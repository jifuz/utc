function [KF_xy]=KF_ai(measurements, model)

%parse input
x=measurements.x;
y=measurements.y;
yaw=measurements.yaw;
dx=measurements.dx;
dy=measurements.dy;
dyaw=measurements.dyaw;
t=measurements.time_s;
t_max=length(x);


%model
Fmodel=[1 0 0 1 0 0;...
        0 1 0 0 1 0;...
        0 0 1 0 0 1;...
        0 0 0 1 0 0;...
        0 0 0 0 1 0;...
        0 0 0 0 0 1];
Q=eye(6).*0.0001;
R=eye(2);



%initial M and V
M=[x(1); y(1); yaw(1); dx(1); dy(1) dyaw(1)];
V=eye(6);
MM=cell(t_max,1);
VV=cell(t_max,1);
MM{1}=M;
VV{1}=V;

%KF main loop
for i=2:t_max
    
   %predict
    [M,V]=ekf_predict(M,V,Fmodel,Q,i);
    
    %update
    [M,V]=ekf_update(M,V,Gmodel,measurements,R,i);
    
    MM{i}=M;
    VV{i}=V;
    
end




end



function [M_,V_]=ekf_predict(M,V,Fmodel,Q,i)

%M
M_=Fmodel*M;

%V
dF=Fmodel;
V_=dF*V*dF'+Q;

end

 



function [M_,V_]=ekf_update(M,V,Gmodel,measurements,R,i)

Gm=driveOnNearestRoadSeg(Gmodel,M);

mx=measurements.x;
my=measurements.y;
measurement=[mx(i);my(i)];
Gm=C*M;

innovation=(measurement-Gm);


Jg=C;


gain=V*Jg' / (Jg*V*Jg'+R);
M_=M+gain*innovation;
V_=(eye(6)-gain*Jg)*V;

end

function Gm=driveOnNearestRoadSeg(Gmodel,M)

lane=Gmodel.lane;


end