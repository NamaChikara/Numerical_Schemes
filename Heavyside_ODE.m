% Note that we are not actually using the Heavyside graph - it is not a
% function at x=0.  Instead, we use {(x,0) for x<0}U{(x,1) for x>=0}.

% Our ODE is u'+A(u)+f(t)=0, so explicit Euler is
%     u(n+1)=u(n)+h*(-A(u(n))-f(t(n)))

h=0.1;  % time step
t=[0:h:10];   % time grid

f=@(t) sin(t)-0.5; 
IC=2;    % initial value

% Explicit Euler
u=zeros(1,length(t)); % initialize u solution vector for explicit Euler
du=zeros(1,length(t)-1); % du is calculated for previous u value
tdu=t(1:length(t)-1);    %    in the loop below. no reason to have
                         %    du(length(t)) value as it wont be used
u(1)=IC;
for i=2:length(t)
    % t(1)=u(1), so we move to t(2),u(2)
    if u(i-1)<0
        du(i-1)=-f(t(i-1));
        u(i)=u(i-1)+h*du(i-1);
    else du(i-1)=-1+-f(t(i-1));   % see note at beginning of code
        u(i)=u(i-1)+h*du(i-1);
    end
end

% Implict Euler with Resolvent
Ru=zeros(1,length(t)); 
Ru(1)=IC;
for i=2:length(t)
    Ru(i)=Resolvent_Hvsd(Ru(i-1),f(t(i)),h);
end

sol=zeros(1,length(t)); % analytical solution
a=zeros(length(t));
for i=1:length(t)
    u1=u(i);
    if i==1
        u0=u(i);
        a(i)=0;
    else u0=u(i-1);
    end
    if sign(u1)-sign(u0)>0
        a(i)=t(i);
    else if sign(u1)-sign(u0)<0
        a(i)=-t(i);
        end
    end
    if u(i)<0
        %sol(i)=log(abs(cos(t(i))))+IC+0+sum(a(1:i));
        sol(i)=cos(t(i))+0.5*t(i)+IC-1+sum(a(1:i));  % -F(t)+u_0+F(0)+c
    else %sol(i)=log(abs(cos(t(i))))-t(i)+IC+0+sum(a(1:i));
        sol(i)=cos(t(i))+0.5*t(i)-t(i)+IC-1+sum(a(1:i)); %-F(t)-t+u(0)+F(0)+c
    end
end

figure
plot(t,u,'r',t,Ru,'-o',t,sol,'--')%tdu,du,'k',
legend('Explicit Euler','Resolvent Method','Analytical')
%title('Solutions for f(t)=sin(t)-0.5 with u(0)=2')
xlabel('t');ylabel('u(t)')
    