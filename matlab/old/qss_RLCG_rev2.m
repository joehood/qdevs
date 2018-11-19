clear all

R = 1;
L = 1;
E = 1;
G = 1;
C = 1;
H = 1;

% State space benchmark
%----------------------------------------
h = 0.1;
A = [-R/L  -1/L ; 1/C  -G/C ];
B = [ 1/L 0 ; 0 1/C];
u = [ E  H ]';

Apr = inv(eye(2) - h * A)
Bpr = h * B * Apr;

tstop = 5.0;

t = 0:h:tstop;
n = length(t);
x = zeros(2, n);

for k=2:n
    x(:,k) = Apr * x(:,k-1) + Bpr * u; 
end

% -----------------------------------------

% history:
thist1(1) = 0.0;
qhist1(1) = 0.0;
k1 = 1;
thist2(1) = 0.0;
qhist2(1) = 0.0;
k2 = 1;

function x = dint(t, tlast, x0, d)
    dt = t - tlast;
    x = x0 + d * dt; 
end

function tnext = ta(d, x, q, dQ, t)
    if (d > 0)
        tnext = t + (q + dQ - x) / d;
    elseif (d < 0)
        tnext = t + (q - 0.5 * dQ - x) / d;
    else
        tnext = inf;
    end
end

function d1 = f1(R, L, E, q1, q2)
    d1 = 1/L * (E - q1*R - q2);
end

function d2 = f2(G, C, H, q1, q2)
    d2 = 1/C * (H - q2*G + q1);
end

function q = quantize(x, q, dQ, eps)

    if x >= q + dQ - eps
        q = q + dQ;
    elseif x <= q - 0.5 * dQ + eps
        q = q - dQ;
    end

end

function [t, tnext1, tlast1, x1, d1, q1, qlast1] = update1(tnext1, tlast1, x1, x2, d1, q1, qlast1, dQ1, eps1, R, L, E)

    global thist1, qhist1, k1
    
    t = tnext1;
    x1 = dint(t, tlast1, x1, d1);
    d1 = f1(R, L, E, x1, x2);
    tnext1 = ta(d1, x1, q1, dQ1, t);
    q1 = quantize(x1, q1, dQ1, eps1);
    tlast1 = t; 
    
    if q1 ~= qlast1
        qlast1 = q1;
        thist1(k1) = t;
        qhist1(k1) = q1;
        k1 = k1 + 1;
        tnext2 = t;
    end
end

function [t, tnext2, tlast2, x2, d2, q2, qlast2] = update2(tnext2, tlast2, x1, x2, d2, q2, qlast2, dQ2, eps2, G, C, H)

    global thist2, qhist2, k2
    
    t = tnext2;
    x2 = dint(t, tlast2, x2, d2);
    d2 = f2(G, C, H, x1, x2);
    tnext2 = ta(d2, x2, q2, dQ2, t);
    q2 = quantize(x2, q2, dQ2, eps2);
    tlast2 = t;
    if q2 ~= qlast2
        qlast2 = q2;
        thist2(k2) = t;
        qhist2(k2) = q2;
        k2 = k2 + 1;
        tnext1 = t;
    end
    
end

% basic QSS solution:

% states for model 1:
x1 = 0.0;      % current 
d1 = 0.0;      % derivative for x1
q1 = 0.0;      % quantized x1
qlast1 = 0.0;  % last quantized x1
dQ1 = 0.01;    % quantum for x1
eps1 = dQ1/4;  % his. width for x1
tlast1 = 0.0;  % previous update time for x1
tnext1 = inf;  % tnext for x1 

% states for model 2:
x2 = 0.0;      % voltage
d2 = 0.0;      % derivative for x2
q2 = 0.0;      % quantized x2
qlast2 = 0.0;  % last quantized x2
dQ2 = 0.01;    % quantum for x2
eps2 = dQ2/4;  % his. width for x2
tlast2 = 0.0;  % previous update time for x2
tnext2 = inf;  % tnext for x2

% initialize simulator:
t = 0.0;

% initialize models:
d1 = f1(R, L, E, q1, q2);
d2 = f2(G, C, H, q1, q2);

% initialize tnext:
tnext1 = ta(d1, x1, q1, dQ1, t);
tnext2 = ta(d2, x2, q2, dQ2, t);

while (t < tstop)
    
    if tnext1 == tnext2
        [t, tnext1, tlast1, x1, d1, q1, qlast1] = update1(tnext1, tlast1, x1, d1, q1, qlast1, dQ1, eps1, R, L, E);
        [t, tnext2, tlast2, x2, d2, q2, qlast2] = update2(tnext2, tlast2, x2, d2, q2, qlast2, dQ2, eps2, G, C, H);
    elseif tnext1 < tnext2       
        [t, tnext1, tlast1, x1, d1, q1, qlast1] = update1(tnext1, tlast1, x1, d1, q1, qlast1, dQ1, eps1, R, L, E);
    elseif tnext2 < tnext1    
        [t, tnext2, tlast2, x2, d2, q2, qlast2] = update2(tnext2, tlast2, x2, d2, q2, qlast2, dQ2, eps2, G, C, H);       
    end
    
end

plot(t, x(1,:), 'c--'); hold on;
plot(thist1, qhist1, 'b-o');
plot(t, x(2,:), 'm--');
plot(thist2, qhist2, 'r-o');

%disp('done.')

