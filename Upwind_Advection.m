function err=Upwind_Advection(a,tpoints)

% Method is only stable if |ak/h|\leq 1
% Recall that h=delta(x), k=delta(t)

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

vUp=a*k/h;  % constant for upwind

if (abs(a*k/h)>1)
	error('Please select a smaller a.')
end

%u1=sin(x);	% u(x,0), given
u1=zeros(1,length(x));
u1(1:2*round(length(x)/xe))=1; % u(x,0)=1 if 0<=x<2
u1(2*round(length(x)/xe)+1:end)=0;  % u(x,0)=0 if 2<x

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
    plot(x,sol_matrix_Up(k,:),'-') 
    pause(0.05);     
end
