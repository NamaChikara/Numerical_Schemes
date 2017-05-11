datapoints=6;
o=2;  % desired rate of convergence

h=zeros(1,datapoints);
E1=zeros(1,datapoints);
E2=zeros(1,datapoints);
p1=zeros(1,datapoints-1);
p2=zeros(1,datapoints-1);

h(1)=1.6667e-1;
h(2)=5.8823e-2;
h(3)=2.1266e-2;
h(4)=1.0753e-2;
h(5)=7.1942e-3;
h(6)=4.7393e-3;

E1(1)=3.5861e-4;
E1(2)=6.4747e-5;
E1(3)=8.1337e-6;
E1(4)=2.3499e-6;
E1(5)=5.5261e-7;
E1(6)=2.4354e-7;

E2(1)=1.1765e-17;
E2(2)=3.6872e-4;
E2(3)=4.2115e-4;
E2(4)=4.2743e-4;
E2(5)=4.2918e-4;
E2(6)=4.2948e-4;

for i=2:datapoints;    
    p1(i)=(log(E1(i))-log(E1(i-1)))/(log(h(i))-log(h(i-1)));
    p2(i)=(log(E2(i))-log(E2(i-1)))/(log(h(i))-log(h(i-1)));
end

p1avg=sum(p1)/length(p1)
p2avg=sum(p2)/length(p2)

loglog(h,E1,'r--o',h,E2,'b--o',h,h.^o,'k')
grid on
xlabel('Mesh Size')
ylabel('Error')
lgn=legend('Lin. Ana. vs Lin. Num.','Ana. vs Lin. Num.','Order 2 Convergence');
lgn.Location='west';
hold off

%{
Sp17_R4_T1  (used Rate_Of_Convergence v5.10.n1)

Rates of Convergence using Two_Point_BVP_Smoothing v5.10.n1
E1 is the 2-norm grid error of Linear Analytical vs. Linear Numerical
E2 is Linear Numerical vs. (unsmoothed) Analytical

h(1)=1.6667e-1;
h(2)=5.8823e-2;
h(3)=2.1266e-2;
h(4)=1.0753e-2;
h(5)=7.1942e-3;
h(6)=4.7393e-3;

E1(1)=3.5861e-4;
E1(2)=6.4747e-5;
E1(3)=8.1337e-6;
E1(4)=2.3499e-6;
E1(5)=5.5261e-7;
E1(6)=2.4354e-7;

E2(1)=1.1765e-17;
E2(2)=3.6872e-4;
E2(3)=4.2115e-4;
E2(4)=4.2743e-4;
E2(5)=4.2918e-4;
E2(6)=4.2948e-4;

%}
