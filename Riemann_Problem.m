function Riemann_Problem(a,tpoints)

if (a==0)
	error('Please select a nonzero a.')
end

% System Parameters

l=90; % number of x grid points
xe=30; % final x value
h=xe/l; % dx
x=[0:h:xe]; % x grid

m=tpoints; % number of t grid points
te=2; % final t value
k=te/m; % dt
t=[0:k:te]; % t grid

if (abs(a*k/h)>1)
	error(['Please select "a" smaller than ', num2str(h/k), '. ',...
        num2str(0.95*h/k), ' gives |ak/h|=0.95.']) % Upwind stable for |ak<h|<1
end

% Exact Solution

sol_mat_exact=zeros(tpoints,length(x));   % Matrix whose i'th row is 
exact_u=zeros(1,length(x));               %     the exact solution at
                                          %     time t(i)
for i=1:length(t)
    for j=1:length(x)
        if j*h<2+a*i*k
            exact_u(j)=1;
        else exact_u(j)=0;
        end
    end
    sol_mat_exact(i,:)=exact_u;
end

% Lax-Friedrichs

vLax=a*k/(2*h); % constant for Lax-Friedrichs method

u1=zeros(1,length(x));
u1(1:2*round(length(x)/xe))=1; % u(x,0)=1 if 0<=x<2
u1(2*round(length(x)/xe)+1:end)=0;  % u(x,0)=0 if 2<x

sol_mat_Lax=zeros(tpoints,length(x));
sol_mat_Lax(1,:)=u1;

for i=2:length(t)
    u0=u1;
    u1(1)=1;
    u1(length(x))=0;
    u1(2:length(x)-1)=(0.5+vLax)*u0(1:length(x)-2)+(0.5-vLax)*u0(3:length(x));
    sol_mat_Lax(i,:)=u1;
end

Err_Lax=sqrt(h).*sqrt(sum(abs(sol_mat_exact-sol_mat_Lax).^2,2));

% Upwind 

vUp=a*k/h;  % constant for upwind

u1=zeros(1,length(x));
u1(1:2*round(length(x)/xe))=1; % u(x,0)=1 if 0<=x<2
u1(2*round(length(x)/xe)+1:end)=0;  % u(x,0)=0 if 2<x

sol_mat_Up=zeros(tpoints,length(x));
sol_mat_Up(1,:)=u1;

for i=2:length(t)
	u0=u1;	% update (i-1)k solution vector
	u1(2:l+1)=u0(2:l+1)-vUp*(u0(2:l+1)-u0(1:l));  % upwind for a>0
	u1(1)=1; 	% Reimann condition
    sol_mat_Up(i,:)=u1;
end

Err_Up=sqrt(h).*sqrt(sum(abs(sol_mat_exact-sol_mat_Up).^2,2));

% Plot Solutions

figure
for k = 1:tpoints
    plot(x,sol_mat_Lax(k,:),'r',[0 xe],[Err_Lax(k) Err_Lax(k)],'r--',...
        x,sol_mat_Up(k,:),'b',[0 xe],[Err_Up(k) Err_Up(k)],'b--',...
        x,sol_mat_exact(k,:),'k');
    pause(0.05);     
end
    title('Riemann Problem: Comparison of Numerical Solutions')
    lgnd=legend('Lax-F','Lax-F error','Upwind','Upwind error','Exact');
    lgnd.Location='northwest';

% Plot Error

%figure
%plot(t,Err_Lax','r--',t,Err_Up,'b--')
