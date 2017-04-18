function v=Yosida_Hvsd(a,x)
% a is Yosida parameter
% x is input variable

if x<0
    v=0;
else if x<a
        v=(1/a)*x;
    else v=1;
    end
end
 