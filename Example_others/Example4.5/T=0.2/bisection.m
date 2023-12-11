function c = bisection(xi,time)
a = -1.0;
b = 2.0;
eps = 1e-15;
c_d = b-a;
while (abs(c_d) > eps)
    c = (a+b)/2.0;
    fa = a+sin(2.0*pi*a)*time-xi;
    fc = c+sin(2.0*pi*c)*time-xi;
    if (fa*fc < 0.0)
        b = c;
        c_d = b-a;
    else
        a = c;
        c_d = b-a;
    end
end
end