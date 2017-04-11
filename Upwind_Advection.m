function error=Upwind_Advection(a)

% Method is only stable if |ak/h|\leq 1
% Recall that h=delta(x), k=delta(t)

if (a==0)
	error('Please select a nonzero a.')
end

m=20;		% # time steps
l=20;		% # x steps
k=0.1;		% dt
h=2*pi/l;	% dx
x=[0:h:2*pi] % l+1 length x-grid

if (abs(a*k/h)>1)
	error('Please select a smaller a.')
end

u1=sin(x)	% u(x,0), given

if a>0
	then for i=2:m
	u0=u1;	% update (i-1)k solution vector
	u1(2:l+1)=u0(2:l+1)-(a*k/h)*(u0(2:l+1)-u0(1:l));  % upwind for a>0
	u1(1)=u1(l+1); 	% periodic BC prescribed
end

else 
	for i=2:m
	u0=u1;
	u1(1:l)=u0(1:l)-(a*k/h)*(u0(2:l+1)-u0(1:l));  % upwind for a<0
	u1(l+1)=u1(1);	% periodic BC
end

plot(x,u1,'-*')
hold on
plot(x,sin(x-a*k*m),'k')  % exact solution: sin(x-a*t_final)
title('Upwind Scheme vs. Exact Solution')

error=norm(sin(x-a*k*m)-u1);
end