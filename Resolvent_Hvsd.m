function out=Resolvent_Hvsd(un,fn1,tau) % u(n), f(t(n+1)), time step/Yosida parameter

x=un-tau*fn1; % what the Resolvent is evaluating

if x<0
    out=x;
else if (0<=x)&&(x<=tau)
        out=0;
    else out=x-tau;
    end
end

