function xyz=leapfrog17v2(a)

% a is the scalar parameter in our scalar advection equation,
%   u_t+au_x=0.
% h=delta(x), k=delta(t)

% Note: |ak/h|<1 for stability 

k=0.1;
tsteps=20;   % number of time grid points
m=19;        % number of grid points - 1
h=2*pi/(m+1);
x=[0:h:2*pi];

% The leapfrog method is given by 
%   U(jh,(n+1)k)=U(jh,(n-1)k)-(ak/h)(U((j+1)h,nk)-U((j-1)h,nk))
% where U(jh,nk)~u(jh,nk)=exact solution at (x,t)=(jh,nk).

% Initial conditions for  u(x,0) and  u(x,k)  are required.
% We can either generate  u(x,k)  from a different method, or
%   cheat and use the exact solution.
% If  u(x,0)=f(x),  u(x,t)=f(x-at).

u1=sin(x);      % u(x,0)
u2=sin(x-a*k);  % u(x,k) from exact solution  (x,t)=sin(x-a*t)

for n=2:tsteps 
u0=u1;  % This is necessary b/c the function below updates the next
u1=u2;  %   time step by using the previous two times.  The most recent
        %   time is the u2 vector from the previous iteration.  If we
        %   rewrite it as we go along, we will lose the grid points we
        %   need.  So  previous time step -> u1,  previous previous time
        %   step -> u2, store values into u2 as current time

% grid points for t=2, x=[2h,(m+1)h]
   u2(2:m+1)=u0(2:m+1)-(a*k/h)*(u1(3:m+2)-u1(1:m));

% periodic boundary conditions, u(0,t)=u((m+2)h,t)
   u2(1)=u0(1)-(a*k/h)*(u1(2)-u1(m+2));
   u2(m+2)=u2(1);
end

% plot result at  t=20*k  vs  exact solution
plot(x,u2,'k',x,sin(x-a*k*tsteps),'b')
xlabel('x')
ylabel('u')
title('Leapfrog Solution vs. Exact Solution')