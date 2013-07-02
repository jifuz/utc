%% compute road boundary
function [xl,yl,xr,yr]=roadBoundary(x,y,w)

xl=x;
yl=y;
xr=x;
yr=y;
numNodes=length(x);

for i=1:numNodes-1
    p0=[x(i) y(i)];
    p1=[x(i+1) y(i+1)];
    
    vertV=vertVector(p1-p0);
    
    p0l=p0+vertV*w/2;
    p0r=p0-vertV*w/2;
    
    xl(i)=p0l(1);
    yl(i)=p0l(2);
    xr(i)=p0r(1);
    yr(i)=p0r(2);
    
end

%last point pair
p0=[x(i) y(i)];
p1=[x(i+1) y(i+1)];
vertV=vertVector(p1-p0);
p1l=p1+vertV*w/2;
p1r=p1-vertV*w/2;

xl(i+1)=p1l(1);
yl(i+1)=p1l(2);
xr(i+1)=p1r(1);
yr(i+1)=p1r(2);


end

function vertG=vertVector(G)
vertG=-[G(2) -G(1)];
vertG=vertG/norm(vertG);
end