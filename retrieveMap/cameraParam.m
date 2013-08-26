function Cam=cameraParam(lat,lon,G)
nx=240;%375;%number of pixel (x)
ny=795;%1242;%number of pixel (y)
n=-1;%near plane distance z
fov=70*2;%fov=30.88*2;%Field-Of-View
f=-100;%far plane distance z

E=[0 0 1.65*6.7];%1.2*15];%Camera center
%G=[10 2 0];%Gaze direction
T=[0 0 1];%View-up vector

Cam.nx=nx;
Cam.ny=ny;
Cam.n=n;
Cam.fov=fov;
Cam.f=f;
Cam.E=E;
Cam.G=G;
Cam.T=T;

Cam.Gl=[rotateVectorClockwise(G(1:2),fov/2); 0]';
Cam.Gr=[rotateVectorClockwise(G(1:2),-fov/2); 0]';

Cam.lat=lat;
Cam.lon=lon;
[latUTM,lonUTM]=deg2utm(lat,lon);
Cam.latUTM=latUTM;
Cam.lonUTM=lonUTM;
end