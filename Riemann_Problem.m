function out=Riemann_Problem(a)

l=50; % number of x grid points
xe=10; % final x value
h=l/xe; % dx
x=[0:h:xe]; % x grid

m=20; % number of t grid points
te=10; % final t value
k=te/m; % dt
t=[0:k:te]; % t grid

v=a*k/(2*h); % value for Lax-Friedrichs method

u=zeros(1,length(x));
u0=zeros(1,length(x));
u0(1:round(length(x)/xe))=1; % u(x,0)=1 if 0<=x<1
u0(round(length(x)/xe)+1:end)=0;  % u(x,0)=0 if 1<x

for i=1:length(t)
    u(1)=1;
    u(length(x))=0;
    for j=2:length(x)-1
        u(j)=(0.5+v)*u0(j-1)+(0.5-v)*u0(j+1);
    end
    u0=u;
end

plot(x,u0)
out=4;