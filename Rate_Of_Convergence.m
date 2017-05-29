datapoints=5;
o=2;  % desired rate of convergence

h=zeros(1,datapoints);
E1=zeros(1,datapoints);
%E2=zeros(1,datapoints);
p1=zeros(1,datapoints-1);
%p2=zeros(1,datapoints-1);

h(1)=10e-1;
h(2)=10e-2;
h(3)=10e-3;
h(4)=10e-4;
h(5)=10e-5;

E1(1)=4.297765067677046e-04;
E1(2)=8.644868461637017e-06;
E1(3)=7.217035286021560e-06;
E1(4)=7.216878437104479e-06;
E1(5)=7.216878437094887e-06;

for i=2:datapoints;    
    p1(i)=(log(E1(i))-log(E1(i-1)))/(log(h(i))-log(h(i-1)));
    %p2(i)=(log(E2(i))-log(E2(i-1)))/(log(h(i))-log(h(i-1)));
end

p1 

loglog(h,E1,'r--o',h,h.^o,'k')  % ...h,E2,'b--o',...
grid on
xlabel('Regularization Width')
ylabel('Error')
lgn=legend('Numerical vs Linear Numerical','Order 2 Convergence');
lgn.Location='west';
hold off

%{
Sp17_R5_T1  (used Rate_Of_Convergence v5.22.n1)

NOTE: h(i) are actually epsilon values in this case. see * below

Rates of Convergence using Two_Point_BVP_Smoothing v5.22.n1
E1 is the 2-norm grid error of numerical solutions to the non-regularized
problem vs the numerical solutions to the linearly regularized problem.

* We hold h constant while refining our intervals of regularization,
epsilon.  The epsilon refinements are h(i)=10e-i.  We refine the grid 5
times,  grid j width = 10e-j.  The h data is included first, then each set
of errors is labeled with the respective grid width we are using.

epsilon refinements:
h(1)=10e-1;
h(2)=10e-2;
h(3)=10e-3;
h(4)=10e-4;
h(5)=10e-5;

error with grid width=10e-1:
E1(1)=0.007288689868557;
E1(2)=0.007288689868557;
E1(3)=0.007288689868557;
E1(4)=0.007288689868557;
E1(5)=0.007288689868549;

error with grid width 10e-2:

E2(1)=8.380824690923966e-04;
E2(2)=7.217600016625979e-04;
E2(3)=7.217600016625990e-04;
E2(4)=7.217600016626800e-04;
E2(5)=7.217600016618783e-04;

error with grid width 10e-3:

E3(1)=4.356973803768319e-04;
E3(2)=7.232256229193578e-05;
E3(3)=7.216885581742335e-05;
E3(4)=7.216885581743114e-05;
E3(5)=7.216885581735129e-05;

error with grid with 10e-4:

E4(1)=4.297765067677046e-04;
E4(2)=8.644868461637017e-06;
E4(3)=7.217035286021560e-06;
E4(4)=7.216878437104479e-06;
E4(5)=7.216878437094887e-06;
%}

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


