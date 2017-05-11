% We're considering the second order ODE with two BCs given by
%   u''(x)= 0 if x<0.5;  1 if x>0.5,  u(0)=u(1)=0.
% We're using the Centered Difference Approximation (Leveque pg. 15-16)
%   D^2U_j=\frac{U_{j-1}-2U_j+U_{j+1}}{h^2}.
% Since u'' is not differentiable, the solution, u, is only twice
%   differentiable.  To do h^2 accuracy analysis of the numerical scheme,
%   we would like u to be four times differentiable.  
% We first linearize the discontinuity at x=0.5, then check quadratic
% smoothing. At each stage, we compare solutions.
l=211; % number of x grid points
xe=1; % final x value
h=xe/l % dx
x=[0:h:xe]; % x grid

% Numerical solution without smoothing

u=zeros(length(x),1);  % Initialize solution vector
u(1)=0;     % Boundary condition, u(x=0)
u(length(x))=0;   % Boundary condition, u(x=1)

F=zeros(length(x)-2,1);   % Column vector for numerical scheme.
F(1)=0; 
F(end)=1;
for i=2:length(x)-3
    if x(i+1)<0.5
        F(i)=0;
    else F(i)=1;
    end
end

A=zeros(length(x)-2,length(x)-2);   % Matrix for numerical scheme.
A(1,1)=-2;  % Columns 1 and m done manually because they do not
A(1,2)=1;   %   include both the lower and upper diagonals.
A(length(x)-2,length(x)-2)=-2;  
A(length(x)-2,length(x)-3)=1;
for i=2:length(x)-3; 
    A(i,i)=-2;  % main diagonal is -2
    A(i,i-1)=1; % upper diagonal is 1
    A(i,i+1)=1; % lower diagonal is 1
end

u(2:length(x)-1)=h^2.*(A^(-1)*F);
u=u';

% Exact solution without smooting

ex=zeros(1,length(x));
for i=1:length(x)
    if x(i)<0.5
        ex(i)=-0.125*x(i);
    else ex(i)=0.5*x(i)^2-0.625*x(i)+0.125;
    end
end

% Numerical solution to linear smoothing

ul=zeros(length(x),1);  % Initialize solution vector
ul(1)=0;     % Boundary condition, u(x=0)
ul(length(x))=0;   % Boundary condition, u(x=1)

ep=0.1;   % epsilon value

Fl=zeros(length(x)-2,1);   % Column vector for numerical scheme.
Fl(1)=0; 
Fl(end)=1;
for i=2:length(x)-3
    if x(i+1)<0.5-ep
        Fl(i)=0;
    else if (0.5-ep<=x(i+1))&&(x(i+1)<=0.5+ep)
            Fl(i)=(x(i+1)-(0.5-ep))/(2*ep);     
    else Fl(i)=1;
        end
    end
end

ul(2:length(x)-1)=h^2.*(A^(-1)*Fl);
ul=ul';

% Exact solution to linear smoothing  (see Spring 2017, Report 4)

a=(0.5-ep);  % to simplify matrix
b=(0.5+ep);
c=-(2*ep)^(-1);

R=zeros(6,6);
R(1,:)=[1,c,0,0,0,0];
R(2,:)=[0,c,1,0,0,0];
R(3,:)=[a,a*c,0,1,c,0];
R(4,:)=[0,b*c,b,0,c,1];
R(5,:)=[0,0,0,1,0,0];
R(6,:)=[0,0,1,0,0,1];

S=zeros(6,1);
S(1)=(-a^2)/(4*ep);
S(2)=(0.5*b^2-a*b)/(2*ep)-b;
S(3)=(-a^3)/(6*ep);
S(4)=((1/6)*b^3-0.5*a*b^2)/(2*ep)-0.5*b^2;
S(5)=0;
S(6)=-0.5;

C=R^(-1)*S; % coefficient matrix

exl=zeros(1,length(x));
for i=1:length(x)
    if x(i)<0.5-ep
        exl(i)=C(1)*x(i)+C(4);
    else if (0.5-ep<=x(i))&&(x(i)<=0.5+ep)
            exl(i)=((1/6)*x(i)^3-(0.5)*a*x(i)^2+C(2)*x(i)+C(5))*(2*ep)^(-1);
        else exl(i)=0.5*x(i)^2+C(3)*x(i)+C(6);
        end
    end
end

% Grid-Norm Errors
Linear_Analytical_vs_Linear_Numerical=sqrt(h)*norm(exl-ul)
Analytical_vs_Linear_Numerical=sqrt(h)*norm(ex-ul)

% Plots
%plot(x,ex,'k',x,exl,'k--',x,u,'b',x,ul,'b--')  % Analytic, Lin Analytic, Lin Numerical
%legend('Analytical','Linear Analytical','Numerical','Linear Numerical')