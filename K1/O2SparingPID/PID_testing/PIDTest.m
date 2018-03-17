y0 = 0;
K1 = 5;
K2 = 5;
K3 = 0.1;
t = linspace(0,100,100000);
y = 0*t;
y(1) = y0;
for i = 2:length(t)
    deriv = 100*rand() ;
    if i>3
        deriv = deriv + -K3*(y(i-1)-y(i-2))/(t(i-1)-t(i-2)) + K2*trapz(t(1:i-1),-y(1:i-1)) + K1*(-y(i-1));
    end
    y(i) = y(i-1) + (t(i)-t(i-1))*deriv;
end
plot(t,y);