% Consider two chemical species, A and B, who combine to form P when the 
% product of their concentrations (C=A*B) reaches a certain value (S).  
% The related system of ODEs is given as  
%    P'={0 if C<S}{(C-1)^2 if C>=S},
%    A'=f'-P',
%    B'=g'-P',
% where f and g represent the source of A and B.

% Notice that P' can also be expressed as 
%    P'=H(C-S)*(C-1)^2,
% where H is the Heavyside function.  We solve the system for this model of
% P' and also for the Yosida approximation of P',
%    P'=H_a(C-S)*(C-1)^2.


f=@(t)0.2*t; % species A source
df=@(t) 0.2; % species A source rate
g=@(t)0.1*t; % species B source
dg=@(t) 0.1; % species B source rate
S=2; % concentration product required for precipitation

h=0.1; % time step
t=[0:h:20]; % time grid

% Non-Yosida
A=zeros(1,length(t)); % Species A solution 
A(1)=f(t(1)); % Initial value
B=zeros(1,length(t)); % Species B solution
B(1)=g(t(1)); % Initial value
P=zeros(1,length(t)); % Precipitation solution
P(1)=0; % Initial value

du=length(t);
for i=2:length(t)
    C=A(i-1)*B(i-1);    % Check product of concentrations
    if C<S
        du(i-1)=0;      % No precipitation if C<S
    else du(i-1)=(C-1)^2;
    end
    A(i)=A(i-1)+h*(df(t(i-1))-du(i-1)); % A(t)=A(t-1)+h*dA(t-1), 
    B(i)=B(i-1)+h*(dg(t(i-1))-du(i-1)); %   dA(t-1)=dSource(t-1)-dPrec(t-1)
    P(i)=P(i-1)+h*du(i-1);          % No additional source term for P
end

YA=zeros(1,length(t)); % Species A solution
YA(1)=f(t(1)); % Initial value
YB=zeros(1,length(t)); % Species B solution
YB(1)=g(t(1)); % Initial value
YP=zeros(1,length(t)); % Precipitation solution
YP(1)=0; % Initial value

Ydu=length(t);
a=0.5;      % Yosida parameter
for i=2:length(t)
    C=YA(i-1)*YB(i-1);
    if C<S
        Ydu(i-1)=0;
    else if (S<=C)&&(C<=S+a)        % Yosida regularization
            Ydu(i-1)=(1/a)*(C-1)^2;
        else Ydu(i-1)=(C-1)^2;
        end
    end
    YA(i)=YA(i-1)+h*(df(t(i-1))-Ydu(i-1));
    YB(i)=YB(i-1)+h*(dg(t(i-1))-Ydu(i-1));
    YP(i)=YP(i-1)+h*Ydu(i-1);
end

figure
plot(t,A,'r',t,B,'b',t,P,'k',t,YA,'r--',t,YB,'b--',t,YP,'k--')
xlabel('time (s)');ylabel('Concentration (mol/L)');
lgnd=legend('$CO_3^{2-}$','$Ca^{2+}$','$CaCO_3$','$\mathrm{Yosida} \ CO_3^{2-}$','$\mathrm{Yosida} \ Ca^{2+}$','$\mathrm{Yosida} \ CaCO_3$');
set(lgnd,'Interpreter','latex')
lgnd.Location='northwest';

%% v4.24.n2
% Plot 1 
% f=0.2*t, g=0.1*t, S=2, a=1

% Plot 2
% a=0.5

% Plot 3
% a=0.1

% Plot 4
% S=1, a=0.1