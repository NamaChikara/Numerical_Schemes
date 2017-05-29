m=511; % number of unknowns
h=1/(m+1); % mesh width

x=[0:h:1];  % x-grid

u=zeros(m+2,1);  % Initialize solution vector
u(1)=0;     % Boundary condition, u(x=0)
u(m+2)=0;   % Boundary condition, u(x=1)

f=@(x) 0;  % u''=f(x)

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
A(1,1)=2;  % Columns 1 and m done manually because they do not
A(1,2)=-1;   %   include both the lower and upper diagonals.
A(m,m)=2;  
A(m,m-1)=-1;
for i=2:m-1; 
    A(i,i)=2;  % main diagonal is -2
    A(i,i-1)=-1; % upper diagonal is 1
    A(i,i+1)=-1; % lower diagonal is 1
end

u(2:m+1)=h^2.*(A^(-1)*F);
u=u';

%exact=@(x) (-1/pi^2)*sin(pi*x); % exact solution

exact1=@(x) 0.125.*x;
region1=[0:0.01:0.5]; 
exact2=@(x) -0.5.*x.^2+0.625.*x-0.125;
region2=[0.5:0.01:1];


plot(x,u,'r',region1,exact1(region1),'-k',region2,exact2(region2),'-k')
%plot(x,abs(u-exact(x)),'r')

exact=zeros(1,length(x));

for i=1:length(x)
    if x(i)<0.5
        exact(i)=0.125*x(i);
    else exact(i)=-0.5*x(i)^2+0.625*x(i)-0.125;
    end
end

grid_error=sqrt(h)*norm(abs(u-exact))


%%
% v4.30.n2

% Plot 1  (Spring 2017, Report 3)
% Log-Log plot of error vs. h for f(x)=sin(pi*x) and grid refinements of
% 0.5 for 2^{-3} through 2^{-6}.
% x=[2^(-6),2^(-5),2^(-4),2^(-3)]; % mesh size
% y=[1.439e-5,5.752e-5,2.306e-4,9.279e-4]; % associated 2-grid-norm errors
% loglog(x,y,'-o')
% title('Errors in Centered Difference Appx. vs. Mesh Size')
% xlabel('mesh size (h)'); ylabel('Error on (0,1)')

% v4.30.n3
% Plot 2  (Spring 2017, Report 3)
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
% Plot 3   (Spring 2017, Report 3)
% Absolute value of pointwise error for f(x)=(0 if x<0.5)U(1 if x>=0.5) 
%   and grids of m=6,7,14,15 with the exact solution:
%   u(x)=-0.125*x if x<0.5, 0.5*x^2-0.625*x+0.125 if x>=0.5.
% Here, (xi,erri)=(grid with m=i, abs(error) with m=i)
% plot(x6,err6,'r--',x7,err7,'r',x14,err14,'b--',x15,err15,'b')
% lgnd=legend('$h=1/7$','$h=1/8$','$h=1/15$','$h=1/16$');
% set(lgnd,'Interpreter','latex')
% lgnd.Location='northeast';

% v5.6.n1
% Plot 4 (Spring 2017, Report 4)
% Grid-2-Norm error for f(x)=(0 if x<0.5)U(1 if x>=0.5) and grids which do
%   and do not include  x=0.5.  
%     h_even={time steps with x=0.5 included}
%     E_even={respective errors in solution}
% Similarly for H_odd, E_odd.
%h_even=[1.25e-1,4.167e-2,3.125e-2,1.825e-2,1.02e-2,7.46e-3];
%E_even=[5.3488e-3,2.3442e-3,1.8687e-3,1.1952e-3,6.925e-4,5.149e-4];
%h_odd=[1.4286e-1,5.263e-2,2.875e-2,1.887e-2,9.901e-3,7.752e-3];
%E_odd=[1.0333e-2,1.5132e-3,4.546e-4,1.998e-4,5.54e-5,3.40e-5];
% Calculate approximate rate of convergence values
%p_even=zeros(1,length(h_even)-1);
%p_odd=zeros(1,length(h_odd)-1);
%for i=2:length(p_even)+1
%   p_even(i-1)=(log(E_even(i))-log(E_even(i-1)))/(log(h_even(i))-log(h_even(i-1)));
%   p_odd(i-1)=(log(E_odd(i))-log(E_odd(i-1)))/(log(h_odd(i))-log(h_odd(i-1)));
%end
% Average p value
%p_even_avg=(1/length(p_even))*sum(p_even);
%p_odd_avg=(1/length(p_odd))*sum(p_odd);
%Log-Log plot of time step vs. error
%loglog(h_even,E_even,'r--o',h_odd,E_odd,'b--o',h_even,h_even.^2,'k')
%grid on
%legend('Mesh w/ x=0.5','Mesh w/out x=0.5','Order 2 Convergence')
%xlabel('Mesh Width'); ylabel('Error')

% v5.6.n2
% Plot 5
% The "odd" system above yielded order 2 convergence as would be hoped for
%   with the centered difference approximation that we used.  "Even"
%   yielded order 0.85 convergence.  We now test the order of convergence
%   for a set of grids which both do and do not include x=0.5.
% h_mix=[h_even(1),h_odd(2),h_even(3),h_odd(4),h_even(5),h_odd(6)];
% E_mix=[E_even(1),E_odd(2),E_even(3),E_odd(4),E_even(5),E_odd(6)];
% p_mix=zeros(1,length(h_mix)-1);
% for i=2:length(p_mix)+1
%     p_mix(i-1)=(log(E_mix(i))-log(E_mix(i-1)))/(log(h_mix(i))-log(h_mix(i-1)));
% end
% p_mix_avg=(1/length(p_mix))*sum(p_mix)
% loglog(h_mix,E_mix,'r--o',h_mix,h_mix.^2,'k--o')
% grid on
% legend('Mixed Mesh','Order 2 Convergence')
% xlabel('Mesh Width'); ylabel('Error')
