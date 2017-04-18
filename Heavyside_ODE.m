
h=0.1;  % time step
t=[0:h:1];   % time grid

g=@(u) u;   % for heavyside (A)
f=@(t) sin(t)-0.5;  % for du/dt=A+f(t)

u=zeros(1,length(t)); % initialize u solution vector
du=zeros(1,length(t)-1); % 
u(1)=-0.25;    % initial value

for i=2:length(t)
    % t(1)=u(1), so we move to t(2),u(2)
    if u(i-1)<0
        du(i-1)=f(t(i-1));
        u(i)=u(i-1)+h*du(i-1);
    else du(i-1)=1+f(t(i-1));
        u(i)=u(i-1)+h*du(i-1);
    end
end
tdu=t(1:length(t)-1);

sol=zeros(1,length(t)); 





figure
plot(t,u,'r',tdu,du,'k',t,sol,'--')
    