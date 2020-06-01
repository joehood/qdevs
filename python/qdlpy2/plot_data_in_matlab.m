d = load('c:\temp\output_Pm_MW.dat')

t = d(:,1);
x = d(:,2);

plot(t, x)
title("Mechanical Power (MW)")

d = load('c:\temp\output_ip_knots.dat')

t = d(:,1);
x = d(:,2);

figure();
plot(t, x)
title("Speed Speed (knots)")