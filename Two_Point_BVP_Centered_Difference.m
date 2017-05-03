m=15; % number of unknowns
h=1/(m+1); % mesh width

x=[0:h:1];  % x-grid

u=zeros(m+2,1);  % Initialize solution vector
u(1)=0;     % Boundary condition, u(x=0)
u(m+2)=0;   % Boundary condition, u(x=1)

% f=@(x) 0;  % u''=f(x)

F=zeros(m,1);   % Column vector for numerical scheme.
F(1)=f(x(2))-u(1)/h^2;  % f(x=h)-f(x=0)/h^2
F(m)=f(x(m+1))-u(m+2)/h^2;  % f(x=m*h)-f(x=1)/h^2
for i=2:m-1
    if x(i+1)<0.5
        F(i)=0;
    else F(i)=1;
    end
end

A=zeros(m,m);   % Matrix for numerical scheme.
A(1,1)=-2;  % Columns 1 and m done manually because they do not
A(1,2)=1;   %   include both the lower and upper diagonals.
A(m,m)=-2;  
A(m,m-1)=1;
for i=2:m-1; 
    A(i,i)=-2;  % main diagonal is -2
    A(i,i-1)=1; % upper diagonal is 1
    A(i,i+1)=1; % lower diagonal is 1
end

u(2:m+1)=h^2.*(A^(-1)*F);
u=u';

%exact=@(x) (-1/pi^2)*sin(pi*x); % exact solution

exact1=@(x) -0.125.*x;
region1=[0:0.01:0.5]; 
exact2=@(x) 0.5.*x.^2-0.625.*x+0.125;
region2=[0.5:0.01:1];


%plot(x,u,'r',region1,exact1(region1),'-k',region2,exact2(region2),'-k')
%plot(x,abs(u-exact(x)),'r')

grid_error=sqrt(h)*norm(abs(u-exact(x)))


%%
% v4.30.n2

% Plot 1
% Log-Log plot of error vs. h for f(x)=sin(pi*x) and grid refinements of
% 0.5 for 2^{-3} through 2^{-6}.
% x=[2^(-6),2^(-5),2^(-4),2^(-3)]; % mesh size
% y=[1.439e-5,5.752e-5,2.306e-4,9.279e-4]; % associated 2-grid-norm errors
% loglog(x,y,'-o')
% title('Errors in Centered Difference Appx. vs. Mesh Size')
% xlabel('mesh size (h)'); ylabel('Error on (0,1)')

% v4.30.n3
% Plot 2
% Plot of solutions for f(x)=(0 if x<0.5)U(1 if x>=0.5) and grids of
%   m=6,7,14,15 plus the exact solution of 
%   u(x)=-0.125*x if x<0.5, 0.5*x^2-0.625*x+0.125 if x>=0.5
% Here, (xi,ui)=(grid with m=i, solution with m=i)
%   region1=x<0.5, exact1=solution for x<0.5
% plot(x6,u6,'r--',x7,u7,'r',x14,u14,'b--',x15,u15,'b',region1,exact1(region1),'k',region2,exact2(region2),'k')
% lgnd=legend('$h=1/7$','$h=1/8$','$h=1/15$','$h=1/16$','exact');
% set(lgnd,'Interpreter','latex')
% lgnd.Location='north';

% v5.3.n1
% Plot 3
% Absolute value of pointwise error for f(x)=(0 if x<0.5)U(1 if x>=0.5) 
%   and grids of m=6,7,14,15 with the exact solution:
%   u(x)=-0.125*x if x<0.5, 0.5*x^2-0.625*x+0.125 if x>=0.5.
% Here, (xi,erri)=(grid with m=i, abs(error) with m=i)
% plot(x6,err6,'r--',x7,err7,'r',x14,err14,'b--',x15,err15,'b')
% lgnd=legend('$h=1/7$','$h=1/8$','$h=1/15$','$h=1/16$');
% set(lgnd,'Interpreter','latex')
% lgnd.Location='northeast';