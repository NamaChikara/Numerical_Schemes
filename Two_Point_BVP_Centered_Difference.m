m=63; % number of unknowns
h=1/(m+1); % mesh width

x=[0:h:1];  % x-grid

u=zeros(m+2,1);  % Initialize solution vector
u(1)=0;     % Boundary condition, u(x=0)
u(m+2)=0;   % Boundary condition, u(x=1)

f=@(x) sin(pi*x);  % u''=f(x)
F=zeros(m,1);   % Column vector for numerical scheme.
F(1)=f(x(2))-u(1)/h^2;  % f(x=h)-f(x=0)/h^2
F(m)=f(x(m+1))-u(m+2)/h^2;  % f(x=m*h)-f(x=1)/h^2
for i=2:m-1
    F(i)=f(x(i+1));
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

exact=@(x) (-1/pi^2)*sin(pi*x); % exact solution
%plot(x,u,x,exact(x))

grid_error=sqrt(h)*norm(u-exact(x));


%%
% v2.30.n1

% Plot 1
% Log-Log plot of error vs. h for f(x)=sin(pi*x) and grid refinements of
% 0.5 for 2^{-3} through 2^{-6}.
% x=[2e-6,2e-5,2e-4,2e-3]; % mesh size
% y=[1.439e-5,5.752e-5,2.306e-4,9.279e-4]; % associated 2-grid-norm errors
% loglog(x,y)
% title('Errors in Centered Difference Appx. vs. Mesh Size')
% xlabel('mesh size (h)'); ylabel('Error on (0,1)')