m=9; % number of unknowns
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

u=h^2.*(inverse(A)*F);
u=u';

plot(x,u)