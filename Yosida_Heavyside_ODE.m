
h=0.5; % dt
t=[0:h:20]; % time grid

f=@(t) sin(t)-0.5+0.1*t;

u=zeros(1,length(t)); % solution vector
du=zeros(1,length(t)-1); % derivative vector
tdu=t(1:length(t)-1); % for ploting du

u(1)=-1; % initial value

for i=2:length(t)
    du(i-1)=f(t(i-1))+Yosida_Hvsd(1,u(i-1));
    u(i)=u(i-1)+h*du(i-1);
end


figure
plot(t,u,tdu,du,'--')