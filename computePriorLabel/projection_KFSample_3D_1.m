function [xp,yp]=projection(Cam,x,y,z,KFpert)

global currKFpert

% x = [-8 -8 8 8];
% y = [1  20  20 2];
% z = [0  0  0 0];
% 
% nx=100;%number of pixel (x)
% ny=100;%number of pixel (y)
% n=-1;%near plane distance z
% fov=60*2;%Field-Of-View
% f=-100;%far plane distance z
% 
% E=[0 0 10];%Camera center
% G=[10 2 0];%Gaze direction
% T=[0 0 1];%View-up vector

nx=Cam.nx;
ny=Cam.ny;
n=Cam.n;
fov=Cam.fov;
f=Cam.f;
E=Cam.E;
G=Cam.G;
T=Cam.T;


t=tan((fov/2)*pi/180)*abs(n);
r=(nx/ny)*t;

l=-r;
b=-t;


w=-G/norm(G);
u=cross(T,w)/norm(cross(T,w));
v=cross(w,u);


Mo1=[nx/2 0 0 (nx-1)/2;...
    0 ny/2 0 (ny-1)/2;...
    0 0 1 0;...
    0 0 0 1];

Mo2=[2/(r-l) 0 0 0;...
    0 2/(t-b) 0 0;...
    0 0 2/(n-f) 0;...
    0 0 0 1];

Mo3=[1 0 0 -(l+r)/2;...
    0 1 0 -(b+t)/2;...
    0 0 1 -(n+f)/2;...
    0 0 0 1];

Mo=Mo1*Mo2*Mo3;

Mv1=[u(1) u(2) u(3) 0;...
    v(1) v(2) v(3) 0;...
    w(1) w(2) w(3) 0;...
    0    0    0    1];

Mv2=[1 0 0 -E(1);...
    0 1 0 -E(2);...
    0 0 1 -E(3);...
    0 0 0 1];

Mv=Mv1*Mv2;


Mp=[n 0 0 0;...
    0 n 0 0;...
    0 0 n+f -f*n;...
    0 0 1 0];


M=Mo*Mp*Mv;

p1s=zeros(4,length(x));
p1v=zeros(4,length(x));


% tic
p=[x';y';z';ones(1,length(x))];
numSample=length(currKFpert);
ps=cell(1,numSample);
for i=1:numSample
    pert=[currKFpert(i,1); currKFpert(i,2); 0; 0];
    pert=repmat(pert,1,size(p,2));
    
    G=Cam.G;
    pertG=[currKFpert(i,3) currKFpert(i,4)];
    pertG=pertG./norm(pertG);
    pert_t=atan2(pertG(2),pertG(1));
    t=atan2(G(2),G(1));
    t=pert_t-t;
    %rotate clockwise theta degree
    R=[cos(t) -sin(t);...
        sin(t) cos(t)];
    x=reshape(x,1,length(x));
    y=reshape(y,1,length(y));
    xyR=R*[x; y];
    xr=xyR(1,:);
    yr=xyR(2,:);
    p=[xr;yr;z';ones(1,length(x))];
    
    ps{i}=p+pert;
    
end
psMat=cell2mat(ps);
plsMat=M*psMat;
plsMat(1,:)=plsMat(1,:)./plsMat(4,:);
plsMat(2,:)=plsMat(2,:)./plsMat(4,:);
plsMat(3,:)=plsMat(3,:)./plsMat(4,:);

xps=mat2cell(plsMat(1,:)', size(p,2)*ones(1,numSample),1);
yps=mat2cell(plsMat(2,:)', size(p,2)*ones(1,numSample),1);

% tt=toc
% 
% tic
for i=1:length(x)
    p=[x(i); y(i); z(i); 1];
    p1s(:,i)=M*p;
    %p1v(:,i)=Mv*p;
end


p1s(1,:)=p1s(1,:)./p1s(4,:);
p1s(2,:)=p1s(2,:)./p1s(4,:);
p1s(3,:)=p1s(3,:)./p1s(4,:);


% figure(2)
% plot(p1s(1,:),p1s(2,:))
% axis([0 nx 0 ny]) 

xp=p1s(1,:)';
yp=p1s(2,:)';
% ttt=toc
% 
% slower=tt/ttt

xp=xps;
yp=yps;
end