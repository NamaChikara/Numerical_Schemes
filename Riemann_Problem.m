function Riemann_Problem(a,tpoints)

if (a==0)
	error('Please select a nonzero a.')
end

l=90; % number of x grid points
xe=30; % final x value
h=xe/l; % dx
x=[0:h:xe]; % x grid

m=tpoints; % number of t grid points
te=2; % final t value
k=te/m; % dt
t=[0:k:te]; % t grid

if (abs(a*k/h)>1)
	error('Please select a smaller a.')
end

vLax=a*k/(2*h); % constant for Lax-Friedrichs method
vUp=a*k/h;  % constant for upwind

u1=zeros(1,length(x));
u1(1:2*round(length(x)/xe))=1; % u(x,0)=1 if 0<=x<2
u1(2*round(length(x)/xe)+1:end)=0;  % u(x,0)=0 if 2<x

% Lax-Friedrichs

sol_matrix_Lax=zeros(tpoints,length(x));
sol_matrix_Lax(1,:)=u1;

for i=2:length(t)
    u0=u1;
    u1(1)=1;
    u1(length(x))=0;
    u1(2:length(x)-1)=(0.5+vLax)*u0(1:length(x)-2)+(0.5-vLax)*u0(3:length(x));
    sol_matrix_Lax(i,:)=u1;
end

% Upwind 

sol_matrix_Up=zeros(tpoints,length(x));
sol_matrix_Up(1,:)=u1;

for i=2:length(t)
	u0=u1;	% update (i-1)k solution vector
	u1(2:l+1)=u0(2:l+1)-vUp*(u0(2:l+1)-u0(1:l));  % upwind for a>0
	u1(1)=1; 	% Reimann condition
    sol_matrix_Up(i,:)=u1;
end

figure


for k = 1:tpoints
    plot(x,sol_matrix_Lax(k,:),'r') 
    plot(x,sol_matrix_Up(k,:),'-') 
    pause(0.05);     
end