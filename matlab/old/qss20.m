clear all

Rs = 1;
L = 0.01;
E = 5;
Rc = 10;
C = 0.01;
% H = 1;

% State space benchmark
%----------------------------------------
h = 0.001;
A=[-Rs/L,0,0,0,0,0,0,0,0,0,-1/L,0,0,0,0,0,0,0,0,0;
    0,-Rs/L,0,0,0,0,0,0,0,0, 1/L,-1/L,0,0,0,0,0,0,0,0;
    0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0,0,0,0,0,0,0;
    0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0,0,0,0,0,0;
    0,0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0,0,0,0,0;
    0,0,0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0,0,0,0;
    0,0,0,0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0,0,0;
    0,0,0,0,0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0,0;
    0,0,0,0,0,0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L,0;
    0,0,0,0,0,0,0,0,0,-Rs/L,0,0,0,0,0,0,0,0,1/L,-1/L;
    1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0,0,0,0,0,0,0;
    0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0,0,0,0,0,0;
    0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0,0,0,0,0;
    0,0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0,0,0,0;
    0,0,0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0,0,0;
    0,0,0,0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0,0;
    0,0,0,0,0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0,0;
    0,0,0,0,0,0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0,0;
    0,0,0,0,0,0,0,0,1/C,-1/C,0,0,0,0,0,0,0,0,-1/(Rc*C),0;
    0,0,0,0,0,0,0,0,0,1/C,0,0,0,0,0,0,0,0,0,-1/(Rc*C)];
  
B = [   1/L ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   ]'; 
u = [ E  ]';

Apr = inv(eye(20) -A*h);
Bpr = Apr* h * B ;

tstop = 15.0;

t = 0:h:tstop;
x = zeros(20, length(t));

for k=2:length(t)
    if t(k)>8
        u=10;
    end
    x(:,k) = Apr * x(:,k-1) + Bpr * u;
    
end
% -----------------------------------------
% figure(1)
% e=eig(A);
% plot(real(e), imag(e), '*r' ,'MarkerSize', 10); hold on;
% figure(2)
% plot(t, x(10,:), 'c--'); hold on;
% plot(t, x(20,:), 'm--');
% basic QSS solution:

% states for model 1:
x1 = 0.0;      % current 
d1 = 0.0;      % derivative for x1
q1 = 0.0;      % quantized x1
qlast1 = 0.0;  % last quantized x1
dQ1 = 0.0001;    % quantum for x1
eps1 = dQ1/4;  % his. width for x1
tlast1 = 0.0;  % previous update time for x1
tnext1 = inf;  % tnext for x1 

% states for model 2:
x2 = 0.0;      % current
d2 = 0.0;      % derivative for x2
q2 = 0.0;      % quantized x2
qlast2 = 0.0;  % last quantized x2
dQ2 = 0.0001;    % quantum for x2
eps2 = dQ2/4;  % his. width for x2
tlast2 = 0.0;  % previous update time for x2
tnext2 = inf;  % tnext for x2

% states for model 3:
x3 = 0.0;      % voltage
d3 = 0.0;      % derivative for x3
q3 = 0.0;      % quantized x3
qlast3 = 0.0;  % last quantized x3
dQ3 = 0.0001;    % quantum for x3
eps3 = dQ3/4;  % his. width for x3
tlast3 = 0.0;  % previous update time for x3
tnext3 = inf;  % tnext for x3

% states for model 4:
x4 = 0.0;      % voltage
d4 = 0.0;      % derivative for x4
q4 = 0.0;      % quantized x4
qlast4 = 0.0;  % last quantized x4
dQ4 = 0.0001;    % quantum for x4
eps4 = dQ4/4;  % his. width for x4
tlast4 = 0.0;  % previous update time for x4
tnext4 = inf;  % tnext for x4

% states for model 5:
x5 = 0.0;      % voltage
d5 = 0.0;      % derivative for x4
q5 = 0.0;      % quantized x4
qlast5 = 0.0;  % last quantized x4
dQ5 = 0.0001;    % quantum for x4
eps5 = dQ5/4;  % his. width for x4
tlast5 = 0.0;  % previous update time for x4
tnext5 = inf;  % tnext for x4

% states for model 6:
x6 = 0.0;      % voltage
d6 = 0.0;      % derivative for x4
q6 = 0.0;      % quantized x4
qlast6 = 0.0;  % last quantized x4
dQ6 = 0.0001;    % quantum for x4
eps6 = dQ6/4;  % his. width for x4
tlast6 = 0.0;  % previous update time for x4
tnext6 = inf;  % tnext for x4

% states for model 7:
x7 = 0.0;      % voltage
d7 = 0.0;      % derivative for x4
q7 = 0.0;      % quantized x4
qlast7 = 0.0;  % last quantized x4
dQ7 = 0.0001;    % quantum for x4
eps7 = dQ7/4;  % his. width for x4
tlast7 = 0.0;  % previous update time for x4
tnext7 = inf;  % tnext for x4

% states for model 8:
x8 = 0.0;      % voltage
d8 = 0.0;      % derivative for x4
q8 = 0.0;      % quantized x4
qlast8 = 0.0;  % last quantized x4
dQ8 = 0.0001;    % quantum for x4
eps8 = dQ8/4;  % his. width for x4
tlast8 = 0.0;  % previous update time for x4
tnext8 = inf;  % tnext for x4

% states for model 9:
x9 = 0.0;      % voltage
d9 = 0.0;      % derivative for x4
q9 = 0.0;      % quantized x4
qlast9 = 0.0;  % last quantized x4
dQ9 = 0.0001;    % quantum for x4
eps9 = dQ9/4;  % his. width for x4
tlast9 = 0.0;  % previous update time for x4
tnext9 = inf;  % tnext for x4

% states for model 10:
x10 = 0.0;      % voltage
d10 = 0.0;      % derivative for x4
q10 = 0.0;      % quantized x4
qlast10 = 0.0;  % last quantized x4
dQ10 = 0.0001;    % quantum for x4
eps10 = dQ10/4;  % his. width for x4
tlast10 = 0.0;  % previous update time for x4
tnext10 = inf;  % tnext for x4

% states for model 11:
x11 = 0.0;      % voltage
d11 = 0.0;      % derivative for x4
q11 = 0.0;      % quantized x4
qlast11 = 0.0;  % last quantized x4
dQ11 = 0.0001;    % quantum for x4
eps11 = dQ11/4;  % his. width for x4
tlast11 = 0.0;  % previous update time for x4
tnext11 = inf;  % tnext for x4

% states for model 12:
x12 = 0.0;      % voltage
d12 = 0.0;      % derivative for x4
q12 = 0.0;      % quantized x4
qlast12 = 0.0;  % last quantized x4
dQ12 = 0.0001;    % quantum for x4
eps12 = dQ12/4;  % his. width for x4
tlast12 = 0.0;  % previous update time for x4
tnext12 = inf;  % tnext for x4

% states for model 13:
x13 = 0.0;      % voltage
d13 = 0.0;      % derivative for x4
q13 = 0.0;      % quantized x4
qlast13 = 0.0;  % last quantized x4
dQ13 = 0.0001;    % quantum for x4
eps13 = dQ13/4;  % his. width for x4
tlast13 = 0.0;  % previous update time for x4
tnext13 = inf;  % tnext for x4

% states for model 14:
x14 = 0.0;      % voltage
d14 = 0.0;      % derivative for x4
q14 = 0.0;      % quantized x4
qlast14 = 0.0;  % last quantized x4
dQ14 = 0.0001;    % quantum for x4
eps14 = dQ14/4;  % his. width for x4
tlast14 = 0.0;  % previous update time for x4
tnext14 = inf;  % tnext for x4

% states for model 15:
x15 = 0.0;      % voltage
d15 = 0.0;      % derivative for x4
q15 = 0.0;      % quantized x4
qlast15 = 0.0;  % last quantized x4
dQ15 = 0.0001;    % quantum for x4
eps15 = dQ15/4;  % his. width for x4
tlast15 = 0.0;  % previous update time for x4
tnext15 = inf;  % tnext for x4

% states for model 16:
x16 = 0.0;      % voltage
d16 = 0.0;      % derivative for x4
q16 = 0.0;      % quantized x4
qlast16 = 0.0;  % last quantized x4
dQ16 = 0.0001;    % quantum for x4
eps16 = dQ16/4;  % his. width for x4
tlast16 = 0.0;  % previous update time for x4
tnext16 = inf;  % tnext for x4

% states for model 17:
x17 = 0.0;      % voltage
d17 = 0.0;      % derivative for x4
q17 = 0.0;      % quantized x4
qlast17 = 0.0;  % last quantized x4
dQ17 = 0.0001;    % quantum for x4
eps17 = dQ17/4;  % his. width for x4
tlast17 = 0.0;  % previous update time for x4
tnext17 = inf;  % tnext for x4

% states for model 18:
x18 = 0.0;      % voltage
d18 = 0.0;      % derivative for x4
q18 = 0.0;      % quantized x4
qlast18 = 0.0;  % last quantized x4
dQ18 = 0.0001;    % quantum for x4
eps18 = dQ18/4;  % his. width for x4
tlast18 = 0.0;  % previous update time for x4
tnext18 = inf;  % tnext for x4

% states for model 19:
x19 = 0.0;      % voltage
d19 = 0.0;      % derivative for x4
q19 = 0.0;      % quantized x4
qlast19 = 0.0;  % last quantized x4
dQ19 = 0.0001;    % quantum for x4
eps19 = dQ19/4;  % his. width for x4
tlast19 = 0.0;  % previous update time for x4
tnext19 = inf;  % tnext for x4

% states for model 20:
x20 = 0.0;      % voltage
d20 = 0.0;      % derivative for x4
q20 = 0.0;      % quantized x4
qlast20 = 0.0;  % last quantized x4
dQ20 = 0.0001;    % quantum for x4
eps20 = dQ20/4;  % his. width for x4
tlast20 = 0.0;  % previous update time for x4
tnext20 = inf;  % tnext for x4


% initialize simulator:
t = 0.0;

% initialize models:
d1 = f1(Rs, L, E, x1, x11);
d2 = f2(Rs, L,x2, x11,x12);
d3 = f3(Rs, L,x3, x12,x13);
d4 = f4(Rs, L,x4, x13,x14);
d5 = f5(Rs, L,x5, x14,x15);
d6 = f6(Rs, L,x6, x15,x16);
d7 = f7(Rs, L,x7, x16,x17);
d8 = f8(Rs, L,x8, x17,x18);
d9 = f9(Rs, L,x9, x18,x19);
d10 = f10(Rs, L,x10, x19,x20);
d11 = f11(Rc, C, x1,x2, x11);
d12 = f12(Rc, C,x2,x3,x12);
d13 = f13(Rc, C,x3,x4,x13);
d14 = f14(Rc, C,x4,x5,x14);
d15 = f15(Rc, C,x5,x6,x15);
d16 = f16(Rc, C,x6,x7,x16);
d17 = f17(Rc, C,x7,x8,x17);
d18 = f18(Rc, C,x8,x9,x18);
d19 = f19(Rc, C,x9,x10,x19);
d20 = f20(Rc, C,x10,x20);

% initialize tnext:
tnext1 = ta(d1, x1, q1, dQ1, t);
tnext2 = ta(d2, x2, q2, dQ2, t);
tnext3 = ta(d3, x3, q3, dQ3, t);
tnext4 = ta(d4, x4, q4, dQ4, t);
tnext5 = ta(d5, x5, q5, dQ5, t);
tnext6 = ta(d6, x6, q6, dQ6, t);
tnext7 = ta(d7, x7, q7, dQ7, t);
tnext8 = ta(d8, x8, q8, dQ8, t);
tnext9 = ta(d9, x9, q9, dQ9, t);
tnext10 = ta(d10, x10, q10, dQ10, t);
tnext11 = ta(d11, x11, q11, dQ11, t);
tnext12 = ta(d12, x12, q12, dQ12, t);
tnext13 = ta(d13, x13, q13, dQ13, t);
tnext14 = ta(d14, x14, q14, dQ14, t);
tnext15 = ta(d15, x15, q15, dQ15, t);
tnext16 = ta(d16, x16, q16, dQ16, t);
tnext17 = ta(d17, x17, q17, dQ17, t);
tnext18 = ta(d18, x18, q18, dQ18, t);
tnext19 = ta(d19, x19, q19, dQ19, t);
tnext20 = ta(d20, x20, q20, dQ20, t);


% history:
thist1(1) = 0.0;
qhist1(1) = 0.0;
k1 = 1;

thist2(1) = 0.0;
qhist2(1) = 0.0;
k2 = 1;

thist3(1) = 0.0;
qhist3(1) = 0.0;
k3 = 1;

thist4(1) = 0.0;
qhist4(1) = 0.0;
k4 = 1;

thist5(1) = 0.0;
qhist5(1) = 0.0;
k5 = 1;

thist6(1) = 0.0;
qhist6(1) = 0.0;
k6 = 1;

thist7(1) = 0.0;
qhist7(1) = 0.0;
k7 = 1;

thist8(1) = 0.0;
qhist8(1) = 0.0;
k8 = 1;

thist9(1) = 0.0;
qhist9(1) = 0.0;
k9 = 1;

thist10(1) = 0.0;
qhist10(1) = 0.0;
k10 = 1;

thist11(1) = 0.0;
qhist11(1) = 0.0;
k11 = 1;

thist12(1) = 0.0;
qhist12(1) = 0.0;
k12 = 1;

thist13(1) = 0.0;
qhist13(1) = 0.0;
k13 = 1;

thist14(1) = 0.0;
qhist14(1) = 0.0;
k14 = 1;

thist15(1) = 0.0;
qhist15(1) = 0.0;
k15 = 1;

thist16(1) = 0.0;
qhist16(1) = 0.0;
k16 = 1;

thist17(1) = 0.0;
qhist17(1) = 0.0;
k17 = 1;

thist18(1) = 0.0;
qhist18(1) = 0.0;
k18 = 1;

thist19(1) = 0.0;
qhist19(1) = 0.0;
k19 = 1;

thist20(1) = 0.0;
qhist20(1) = 0.0;
k20 = 1;

while (t < tstop)
    
    if (t>8)
        E=10;
    end
    
    % internal update function
    if (tnext1 < tnext2) && (tnext1 < tnext3) &&(tnext1 < tnext4) && (tnext1 < tnext5) &&(tnext1 < tnext6) &&(tnext1 < tnext7) &&(tnext1 < tnext8)&& (tnext1 < tnext9) &&(tnext1 < tnext10) && (tnext1 < tnext11) &&(tnext1 < tnext12) && (tnext1 < tnext13) &&(tnext1 < tnext14) && (tnext1 < tnext15) &&(tnext1 < tnext16) && (tnext1 < tnext17) &&(tnext1 < tnext18) && (tnext1 < tnext19) &&(tnext1 < tnext20)
        
%         update1();
        
        t = tnext1;        
        q1 = quantize(x1, q1, dQ1, eps1);
        x1 = dint(t, tlast1, x1, d1);
        d1 = f1(Rs, L, E, x1, x11);      
        tnext1 = ta(d1, x1, q1, dQ1, t);
        tlast1 = t;
        if q1 ~= qlast1
            qlast1 = q1;
            thist1(k1) = t;
            qhist1(k1) = q1;
            k1 = k1 + 1;
            x11 = dint(t, tlast11, x11, d11);
            tlast11 = t;
            d11 = f11(Rc, C, x1,x2, x11);
            tnext11 = ta(d11, x11, q11, dQ11, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
        end
        
        
        
    elseif (tnext2 < tnext1) && (tnext2 < tnext3) &&(tnext2 < tnext4) && (tnext2 < tnext5) &&(tnext2 < tnext6) &&(tnext2 < tnext7) &&(tnext2 < tnext8)&& (tnext2 < tnext9) &&(tnext2 < tnext10) && (tnext2 < tnext11) &&(tnext2 < tnext12) && (tnext2 < tnext13) &&(tnext2 < tnext14) && (tnext2 < tnext15) &&(tnext2 < tnext16) && (tnext2 < tnext17) &&(tnext2 < tnext18) && (tnext2 < tnext19) &&(tnext2 < tnext20)
        
        t = tnext2;        
        q2 = quantize(x2, q2, dQ2, eps2);
        x2 = dint(t, tlast2, x2, d2);
        d2 = f2(Rs, L,x2, x11,x12);
        tnext2 = ta(d2, x2, q2, dQ2, t);
        
        tlast2 = t;
        if q2 ~= qlast2
            qlast2 = q2;
            thist2(k2) = t;
            qhist2(k2) = q2;
            k2 = k2 + 1;
            x11 = dint(t, tlast11, x11, d11);
            tlast11 = t;
            d11 = f11(Rc, C, x1,x2, x11);
            tnext11 = ta(d11, x11, q11, dQ11, t);
            x12 = dint(t, tlast12, x12, d12);
            tlast12 = t;
            d12 = f12(Rc, C,x2,x3,x12);
            tnext12 = ta(d12, x12, q12, dQ12, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext3 < tnext1) && (tnext3 < tnext2) &&(tnext3 < tnext4) && (tnext3 < tnext5) &&(tnext3 < tnext6) &&(tnext3 < tnext7) &&(tnext3 < tnext8)&& (tnext3 < tnext9) &&(tnext3 < tnext10) && (tnext3 < tnext11) &&(tnext3 < tnext12) && (tnext3 < tnext13) &&(tnext3 < tnext14) && (tnext3 < tnext15) &&(tnext3 < tnext16) && (tnext3 < tnext17) &&(tnext3 < tnext18) && (tnext3 < tnext19) &&(tnext3 < tnext20)
       
        t = tnext3;        
        q3 = quantize(x3, q3, dQ3, eps3);
        x3 = dint(t, tlast3, x3, d3);
        d3 = f3(Rs, L,x3, x12,x13);
        tnext3 = ta(d3, x3, q3, dQ3, t);
        
        tlast3 = t;
        if q3 ~= qlast3
            qlast3 = q3;
            thist3(k3) = t;
            qhist3(k3) = q3;
            k3 = k3 + 1;
            x12 = dint(t, tlast12, x12, d12);
            tlast12 = t;
            d12 = f12(Rc, C,x2,x3,x12);
            tnext12 = ta(d12, x12, q12, dQ12, t);
            x13 = dint(t, tlast13, x13, d13);
            tlast13 = t;
            d13 = f13(Rc, C,x3,x4,x13);
            tnext13 = ta(d13, x13, q13, dQ13, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        
        elseif (tnext4 < tnext1) && (tnext4 < tnext2) &&(tnext4 < tnext3) && (tnext4 < tnext5) &&(tnext4 < tnext6) &&(tnext4 < tnext7) &&(tnext4 < tnext8)&& (tnext4 < tnext9) &&(tnext4 < tnext10) && (tnext4 < tnext11) &&(tnext4 < tnext12) && (tnext4 < tnext13) &&(tnext4 < tnext14) && (tnext4 < tnext15) &&(tnext4 < tnext16) && (tnext4 < tnext17) &&(tnext4 < tnext18) && (tnext4 < tnext19) &&(tnext4 < tnext20)
        
        t = tnext4;        
        q4 = quantize(x4, q4, dQ4, eps4);
        x4 = dint(t, tlast4, x4, d4);
        d4 = f4(Rs, L,x4, x13,x14);
        tnext4 = ta(d4, x4, q4, dQ4, t);
        
        tlast4 = t;
        if q4 ~= qlast4
            qlast4 = q4;
            thist4(k4) = t;
            qhist4(k4) = q4;
            k4 = k4 + 1;
            x13 = dint(t, tlast13, x13, d13);
            tlast13 = t;
            d13 = f13(Rc, C,x3,x4,x13);
            tnext13 = ta(d13, x13, q13, dQ13, t);
            x14 = dint(t, tlast14, x14, d14);
            tlast14 = t;
            d14 = f14(Rc, C,x4,x5,x14);
            tnext14 = ta(d14, x14, q14, dQ14, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         
        elseif (tnext5 < tnext1) && (tnext5 < tnext2) &&(tnext5 < tnext3) && (tnext5 < tnext4) &&(tnext5 < tnext6) &&(tnext5 < tnext7) &&(tnext5 < tnext8)&& (tnext5 < tnext9) &&(tnext5 < tnext10) && (tnext5 < tnext11) &&(tnext5 < tnext12) && (tnext5 < tnext13) &&(tnext5 < tnext14) && (tnext5 < tnext15) &&(tnext5 < tnext16) && (tnext5 < tnext17) &&(tnext5 < tnext18) && (tnext5 < tnext19) &&(tnext5 < tnext20)
        
        t = tnext5;        
        q5 = quantize(x5, q5, dQ5, eps5);
        x5 = dint(t, tlast5, x5, d5);
        d5 = f5(Rs, L,x5, x14,x15);
        tnext5 = ta(d5, x5, q5, dQ5, t);
        
        tlast5 = t;
        if q5 ~= qlast5
            qlast5 = q5;
            thist5(k5) = t;
            qhist5(k5) = q5;
            k5 = k5 + 1;
            x14 = dint(t, tlast14, x14, d14);
            tlast14 = t;
            d14 = f14(Rc, C,x4,x5,x14);
            tnext14 = ta(d14, x14, q14, dQ14, t);
            
            x15 = dint(t, tlast15, x15, d15);
            tlast15 = t;
            d15 = f15(Rc, C,x5,x6,x15);
            tnext15 = ta(d15, x15, q15, dQ15, t);
            
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext6 < tnext1) && (tnext6 < tnext2) &&(tnext6 < tnext3) && (tnext6 < tnext4) &&(tnext6 < tnext5) &&(tnext6 < tnext7) &&(tnext6 < tnext8)&& (tnext6 < tnext9) &&(tnext6 < tnext10) && (tnext6 < tnext11) &&(tnext6 < tnext12) && (tnext6 < tnext13) &&(tnext6 < tnext14) && (tnext6 < tnext15) &&(tnext6 < tnext16) && (tnext6 < tnext17) &&(tnext6 < tnext18) && (tnext6 < tnext19) &&(tnext6 < tnext20)
        
        t = tnext6;        
        q6 = quantize(x6, q6, dQ6, eps6);
        x6 = dint(t, tlast6, x6, d6);
        d6 = f6(Rs, L,x6, x15,x16);
        tnext6 = ta(d6, x6, q6, dQ6, t);
        
        tlast6 = t;
        if q6 ~= qlast6
            qlast6 = q6;
            thist6(k6) = t;
            qhist6(k6) = q6;
            k6 = k6 + 1;          
            x15 = dint(t, tlast15, x15, d15);
            tlast15 = t;
            d15 = f15(Rc, C,x5,x6,x15);
            tnext15 = ta(d15, x15, q15, dQ15, t);
            x16 = dint(t, tlast16, x16, d16);
            tlast16 = t;
            d16 = f16(Rc, C,x6,x7,x16);
            tnext16 = ta(d16, x16, q16, dQ16, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        
        elseif (tnext7 < tnext1) && (tnext7 < tnext2) &&(tnext7 < tnext3) && (tnext7 < tnext4) &&(tnext7 < tnext5) &&(tnext7 < tnext6) &&(tnext7 < tnext8)&& (tnext7 < tnext9) &&(tnext7 < tnext10) && (tnext7 < tnext11) &&(tnext7 < tnext12) && (tnext7 < tnext13) &&(tnext7 < tnext14) && (tnext7 < tnext15) &&(tnext7 < tnext16) && (tnext7 < tnext17) &&(tnext7 < tnext18) && (tnext7 < tnext19) &&(tnext7 < tnext20)
        
        t = tnext7;        
        q7 = quantize(x7, q7, dQ7, eps7);
        x7 = dint(t, tlast7, x7, d7);
        d7 = f7(Rs, L,x7, x16,x17);
        tnext7 = ta(d7, x7, q7, dQ7, t);
        
        tlast7 = t;
        if q7 ~= qlast7
            qlast7 = q7;
            thist7(k7) = t;
            qhist7(k7) = q7;
            k7 = k7 + 1;          
            x16 = dint(t, tlast16, x16, d16);
            tlast16 = t;
            d16 = f16(Rc, C,x6,x7,x16);
            tnext16 = ta(d16, x16, q16, dQ16, t);
            x17 = dint(t, tlast17, x17, d17);
            tlast17 = t;
            d17 = f17(Rc, C,x7,x8,x17);
            tnext17 = ta(d17, x17, q17, dQ17, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext8 < tnext1) && (tnext8 < tnext2) &&(tnext8 < tnext3) && (tnext8 < tnext4) &&(tnext8 < tnext5) &&(tnext8 < tnext6) &&(tnext8 < tnext7)&& (tnext8 < tnext9) &&(tnext8 < tnext10) && (tnext8 < tnext11) &&(tnext8 < tnext12) && (tnext8 < tnext13) &&(tnext8 < tnext14) && (tnext8 < tnext15) &&(tnext8 < tnext16) && (tnext8 < tnext17) &&(tnext8 < tnext18) && (tnext8 < tnext19) &&(tnext8 < tnext20)
        
        t = tnext8;        
        q8 = quantize(x8, q8, dQ8, eps8);
        x8 = dint(t, tlast8, x8, d8);
        d8 = f8(Rs, L,x8, x17,x18);
        tnext8 = ta(d8, x8, q8, dQ8, t);
        
        tlast8 = t;
        if q8 ~= qlast8
            qlast8 = q8;
            thist8(k8) = t;
            qhist8(k8) = q8;
            k8 = k8 + 1;          
            x17 = dint(t, tlast17, x17, d17);
            tlast17 = t;
            d17 = f17(Rc, C,x7,x8,x17);
            tnext17 = ta(d17, x17, q17, dQ17, t);
            x18 = dint(t, tlast18, x18, d18);
            tlast18 = t;
            d18 = f18(Rc, C,x8,x9,x18);
            tnext18 = ta(d18, x18, q18, dQ18, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext9 < tnext1) && (tnext9 < tnext2) &&(tnext9 < tnext3) && (tnext9 < tnext4) &&(tnext9 < tnext5) &&(tnext9 < tnext6) &&(tnext9 < tnext7)&& (tnext9 < tnext8) &&(tnext9 < tnext10) && (tnext9 < tnext11) &&(tnext9 < tnext12) && (tnext9 < tnext13) &&(tnext9 < tnext14) && (tnext9 < tnext15) &&(tnext9 < tnext16) && (tnext9 < tnext17) &&(tnext9 < tnext18) && (tnext9 < tnext19) &&(tnext9 < tnext20)
        
        t = tnext9;        
        q9 = quantize(x9, q9, dQ9, eps9);
        x9 = dint(t, tlast9, x9, d9);
        d9 = f9(Rs, L,x9, x18,x19);
        tnext9 = ta(d9, x9, q9, dQ9, t);
        
        tlast9 = t;
        if q9 ~= qlast9
            qlast9 = q9;
            thist9(k9) = t;
            qhist9(k9) = q9;
            k9 = k9 + 1;          
            x18 = dint(t, tlast18, x18, d18);
            tlast18 = t;
            d18 = f18(Rc, C,x8,x9,x18);
            tnext18 = ta(d18, x18, q18, dQ18, t);
            x19 = dint(t, tlast19, x19, d19);
            tlast19 = t;
            d19 = f19(Rc, C,x9,x10,x19);
            tnext19 = ta(d19, x19, q19, dQ19, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext10 < tnext1) && (tnext10 < tnext2) &&(tnext10 < tnext3) && (tnext10 < tnext4) &&(tnext10 < tnext5) &&(tnext10 < tnext6) &&(tnext10 < tnext7)&& (tnext10 < tnext8) &&(tnext10 < tnext9) && (tnext10 < tnext11) &&(tnext10 < tnext12) && (tnext10 < tnext13) &&(tnext10 < tnext14) && (tnext10 < tnext15) &&(tnext10 < tnext16) && (tnext10 < tnext17) &&(tnext10 < tnext18) && (tnext10 < tnext19) &&(tnext10 < tnext20)
        
        t = tnext10;        
        q10 = quantize(x10, q10, dQ10, eps10);
        x10 = dint(t, tlast10, x10, d10);
        d10 = f10(Rs, L,x10, x19,x20);
        tnext10 = ta(d10, x10, q10, dQ10, t);
        
        tlast10 = t;
        if q10 ~= qlast10
            qlast10 = q10;
            thist10(k10) = t;
            qhist10(k10) = q10;
            k10 = k10 + 1;          
            x19 = dint(t, tlast19, x19, d19);
            tlast19 = t;
            d19 = f19(Rc, C,x9,x10,x19);
            tnext19 = ta(d19, x19, q19, dQ19, t);
            x20 = dint(t, tlast20, x20, d20);
            tlast20 = t;
            d20 = f20(Rc, C,x10,x20);
            tnext20 = ta(d20, x20, q20, dQ20, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
             
        
         elseif (tnext11 < tnext1) && (tnext11 < tnext2) && (tnext11 < tnext3)&& (tnext11 < tnext4) && (tnext11 < tnext5) &&(tnext11 < tnext6) &&(tnext11 < tnext7) &&(tnext11 < tnext8)&& (tnext11 < tnext9) &&(tnext11 < tnext10) &&(tnext11 < tnext12) && (tnext11 < tnext13) &&(tnext11 < tnext14) && (tnext11 < tnext15) &&(tnext11 < tnext16) && (tnext11 < tnext17) &&(tnext11 < tnext18) && (tnext11 < tnext19) &&(tnext11 < tnext20)
        
        t = tnext11;        
        q11 = quantize(x11, q11, dQ11, eps11);
        x11 = dint(t, tlast11, x11, d11);
        d11 = f11(Rc, C, x1,x2, x11);
        tnext11 = ta(d11, x11, q11, dQ11, t);
       
        tlast11 = t;
        if q11 ~= qlast11
            qlast11 = q11;
            thist11(k11) = t;
            qhist11(k11) = q11;
            k11 = k11 + 1;
            x1 = dint(t, tlast1, x1, d1);
            tlast1 = t;
            d1 = f1(Rs, L, E, x1, x11);
            tnext1 = ta(d1, x1, q1, dQ1, t);
            x2 = dint(t, tlast2, x2, d2);
            tlast2 = t;
            d2 = f2(Rs, L,x2, x11,x12);
            tnext2 = ta(d2, x2, q2, dQ2, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext12 < tnext1) && (tnext12 < tnext2) &&(tnext12 < tnext3) && (tnext12 < tnext4)&& (tnext12 < tnext5) &&(tnext12 < tnext6) &&(tnext12 < tnext7) &&(tnext12 < tnext8)&& (tnext12 < tnext9) &&(tnext12 < tnext10) && (tnext12 < tnext11) && (tnext12 < tnext13) &&(tnext12 < tnext14) && (tnext12 < tnext15) &&(tnext12 < tnext16) && (tnext12 < tnext17) &&(tnext12 < tnext18) && (tnext12 < tnext19) &&(tnext12 < tnext20)
        
        t = tnext12;        
        q12 = quantize(x12, q12, dQ12, eps12);
        x12 = dint(t, tlast12, x12, d12);
        d12 = f12(Rc, C,x2,x3,x12);
        tnext12 = ta(d12, x12, q12, dQ12, t);
        
        
        tlast12 = t;
        if q12 ~= qlast12
            qlast12 = q12;
            thist12(k12) = t;
            qhist12(k12) = q12;
            k12 = k12 + 1;
            x2 = dint(t, tlast2, x2, d2);
            tlast2 =t;
            d2 = f2(Rs, L,x2, x11,x12);
            tnext2 = ta(d2, x2, q2, dQ2, t);
            x3 = dint(t, tlast3, x3, d3);
            tlast3 =t;
            d3 = f3(Rs, L,x3, x12,x13);
            tnext3 = ta(d3, x3, q3, dQ3, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext13 < tnext1) && (tnext13 < tnext2) &&(tnext13 < tnext3) && (tnext13 < tnext4)&& (tnext13 < tnext5) &&(tnext13 < tnext6) &&(tnext13 < tnext7) &&(tnext13 < tnext8)&& (tnext13 < tnext9) &&(tnext13 < tnext10) && (tnext13 < tnext11)&& (tnext13 < tnext12) &&(tnext13 < tnext14) && (tnext13 < tnext15) &&(tnext13 < tnext16) && (tnext13 < tnext17) &&(tnext13 < tnext18) && (tnext13 < tnext19) &&(tnext13 < tnext20)
        
        t = tnext13;        
        q13 = quantize(x13, q13, dQ13, eps13);
        x13 = dint(t, tlast13, x13, d13);
        d13 = f13(Rc, C,x3,x4,x13);
        tnext13 = ta(d13, x13, q13, dQ13, t);
        
        
        tlast13 = t;
        if q13 ~= qlast13
            qlast13 = q13;
            thist13(k13) = t;
            qhist13(k13) = q13;
            k13 = k13 + 1;
            x3 = dint(t, tlast3, x3, d3);
            tlast3 =t;
            d3 = f3(Rs, L,x3, x12,x13);
            tnext3 = ta(d3, x3, q3, dQ3, t);
            x4 = dint(t, tlast4, x4, d4);
            tlast4 =t;
            d4 = f4(Rs, L,x4, x13,x14);
            tnext4 = ta(d4, x4, q4, dQ4, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext14 < tnext1) && (tnext14 < tnext2) &&(tnext14 < tnext3) && (tnext14 < tnext4)&& (tnext14 < tnext5) &&(tnext14 < tnext6) &&(tnext14 < tnext7) &&(tnext14 < tnext8)&& (tnext14 < tnext9) &&(tnext14 < tnext10) && (tnext14 < tnext11) &&(tnext14 < tnext12) &&(tnext14 < tnext13) && (tnext14 < tnext15) &&(tnext14 < tnext16) && (tnext14 < tnext17) &&(tnext14 < tnext18) && (tnext14 < tnext19) &&(tnext14 < tnext20)
        
        t = tnext14;        
        q14 = quantize(x14, q14, dQ14, eps14);
        x14 = dint(t, tlast14, x14, d14);
        d14 = f14(Rc, C,x4,x5,x14);
        tnext14 = ta(d14, x14, q14, dQ14, t);
        
        
        tlast14 = t;
        if q14 ~= qlast14
            qlast14 = q14;
            thist14(k14) = t;
            qhist14(k14) = q14;
            k14 = k14 + 1;
            x4 = dint(t, tlast4, x4, d4);
            tlast4 =t;
            d4 = f4(Rs, L,x4, x13,x14);
            tnext4 = ta(d4, x4, q4, dQ4, t);
            x5 = dint(t, tlast5, x5, d5);
            tlast5 =t;
            d5 = f5(Rs, L,x5, x14,x15);
            tnext5 = ta(d5, x5, q5, dQ5, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext15 < tnext1) && (tnext15 < tnext2) &&(tnext15 < tnext3) && (tnext15 < tnext4)&& (tnext15 < tnext5) &&(tnext15 < tnext6) &&(tnext15 < tnext7) &&(tnext15 < tnext8)&& (tnext15 < tnext9) &&(tnext15 < tnext10) && (tnext15 < tnext11) && (tnext15 < tnext12) && (tnext15 < tnext13) &&(tnext15 < tnext14)  &&(tnext15 < tnext16) && (tnext15 < tnext17) &&(tnext15 < tnext18) && (tnext15 < tnext19) &&(tnext15 < tnext20)
        
        t = tnext15;        
        q15 = quantize(x15, q15, dQ15, eps15);
        x15 = dint(t, tlast15, x15, d15);
        d15 = f15(Rc, C,x5,x6,x15);
        tnext15 = ta(d15, x15, q15, dQ15, t);
        
        
        tlast15 = t;
        if q15 ~= qlast15
            qlast15 = q15;
            thist15(k15) = t;
            qhist15(k15) = q15;
            k15 = k15 + 1;
            x5 = dint(t, tlast5, x5, d5);
            tlast5 =t;
            d5 = f5(Rs, L,x5, x14,x15);
            tnext5 = ta(d5, x5, q5, dQ5, t);
            x6 = dint(t, tlast6, x6, d6);
            tlast6 =t;
            d6 = f6(Rs, L,x6, x15,x16);
            tnext6 = ta(d6, x6, q6, dQ6, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext16 < tnext1) && (tnext16 < tnext2) &&(tnext16 < tnext3) && (tnext16 < tnext4)&& (tnext16 < tnext5) &&(tnext16 < tnext6) &&(tnext16 < tnext7) &&(tnext16 < tnext8)&& (tnext16 < tnext9) &&(tnext16 < tnext10) && (tnext16 < tnext11) && (tnext16 < tnext12) && (tnext16 < tnext13) &&(tnext16 < tnext14)  &&(tnext16 < tnext15) && (tnext16 < tnext17) &&(tnext16 < tnext18) && (tnext16 < tnext19) &&(tnext16 < tnext20)
        
        t = tnext16;        
        q16 = quantize(x16, q16, dQ16, eps16);
        x16 = dint(t, tlast16, x16, d16);
        d16 = f16(Rc, C,x6,x7,x16);
        tnext16 = ta(d16, x16, q16, dQ16, t);
        
        
        tlast16 = t;
        if q16 ~= qlast16
            qlast16 = q16;
            thist16(k16) = t;
            qhist16(k16) = q16;
            k16 = k16 + 1;
            x6 = dint(t, tlast6, x6, d6);
            tlast6 =t;
            d6 = f6(Rs, L,x6, x15,x16);
            tnext6 = ta(d6, x6, q6, dQ6, t);
            x7 = dint(t, tlast7, x7, d7);
            tlast7 =t;
            d7 = f7(Rs, L,x7, x16,x17);
            tnext7 = ta(d7, x7, q7, dQ7, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext17 < tnext1) && (tnext17 < tnext2) &&(tnext17 < tnext3) && (tnext17 < tnext4)&& (tnext17 < tnext5) &&(tnext17 < tnext6) &&(tnext17 < tnext7) &&(tnext17 < tnext8)&& (tnext17 < tnext9) &&(tnext17 < tnext10) && (tnext17 < tnext11) && (tnext17 < tnext12) && (tnext17 < tnext13) &&(tnext17 < tnext14)  &&(tnext17 < tnext15) && (tnext17 < tnext16) &&(tnext17 < tnext18) && (tnext17 < tnext19) &&(tnext17 < tnext20)
        
        t = tnext17;        
        q17 = quantize(x17, q17, dQ17, eps17);
        x17 = dint(t, tlast17, x17, d17);
        d17 = f17(Rc, C,x7,x8,x17);
        tnext17 = ta(d17, x17, q17, dQ17, t);
        
        
        tlast17 = t;
        if q17 ~= qlast17
            qlast17 = q17;
            thist17(k17) = t;
            qhist17(k17) = q17;
            k17 = k17 + 1;
            x7 = dint(t, tlast7, x7, d7);
            tlast7 =t;
            d7 = f7(Rs, L,x7, x16,x17);
            tnext7 = ta(d7, x7, q7, dQ7, t);
            x8 = dint(t, tlast8, x8, d8);
            tlast8 =t;
            d8 = f8(Rs, L,x8, x17,x18);
            tnext8 = ta(d8, x8, q8, dQ8, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext18 < tnext1) && (tnext18 < tnext2) &&(tnext18 < tnext3) && (tnext18 < tnext4)&& (tnext18 < tnext5) &&(tnext18 < tnext6) &&(tnext18 < tnext7) &&(tnext18 < tnext8)&& (tnext18 < tnext9) &&(tnext18 < tnext10) && (tnext18 < tnext11) && (tnext18 < tnext12) && (tnext18 < tnext13) &&(tnext18 < tnext14)  &&(tnext18 < tnext15) && (tnext18 < tnext16) &&(tnext18 < tnext17) && (tnext18 < tnext19) &&(tnext18 < tnext20)
        
        t = tnext18;        
        q18 = quantize(x18, q18, dQ18, eps18);
        x18 = dint(t, tlast18, x18, d18);
        d18 = f18(Rc, C,x8,x9,x18);
        tnext18 = ta(d18, x18, q18, dQ18, t);
        
        
        tlast18 = t;
        if q18 ~= qlast18
            qlast18 = q18;
            thist18(k18) = t;
            qhist18(k18) = q18;
            k18 = k18 + 1;
            x8 = dint(t, tlast8, x8, d8);
            tlast8 =t;
            d8 = f8(Rs, L,x8, x17,x18);
            tnext8 = ta(d8, x8, q8, dQ8, t);
            x9 = dint(t, tlast9, x9, d9);
            tlast9 =t;
            d9 = f9(Rs, L,x9, x18,x19);
            tnext9 = ta(d9, x9, q9, dQ9, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
         elseif (tnext19 < tnext1) && (tnext19 < tnext2) &&(tnext19 < tnext3) && (tnext19 < tnext4)&& (tnext19 < tnext5) &&(tnext19 < tnext6) &&(tnext19 < tnext7) &&(tnext19 < tnext8)&& (tnext19 < tnext9) &&(tnext19 < tnext10) && (tnext19 < tnext11) && (tnext19 < tnext12) && (tnext19 < tnext13) &&(tnext19 < tnext14)  &&(tnext19 < tnext15) && (tnext19 < tnext16) &&(tnext19 < tnext17) && (tnext19 < tnext18) &&(tnext19 < tnext20)
        
        t = tnext19;        
        q19 = quantize(x19, q19, dQ19, eps19);
        x19 = dint(t, tlast19, x19, d19);
        d19 = f19(Rc, C,x9,x10,x19);
        tnext19 = ta(d19, x19, q19, dQ19, t);
        
        
        tlast19 = t;
        if q19 ~= qlast19
            qlast19 = q19;
            thist19(k19) = t;
            qhist19(k19) = q19;
            k19 = k19 + 1;
            x9 = dint(t, tlast9, x9, d9);
            tlast9 =t;
            d9 = f9(Rs, L,x9, x18,x19);
            tnext9 = ta(d9, x9, q9, dQ9, t);
            x10 = dint(t, tlast10, x10, d10);
            tlast10 =t;
            d10 = f10(Rs, L,x10, x19,x20);
            tnext10 = ta(d10, x10, q10, dQ10, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
        
        elseif (tnext20 < tnext1) && (tnext20 < tnext2) &&(tnext20 < tnext3) && (tnext20 < tnext4)&& (tnext20 < tnext5) &&(tnext20 < tnext6) &&(tnext20 < tnext7) &&(tnext20 < tnext8)&& (tnext20 < tnext9) &&(tnext20 < tnext10) && (tnext20 < tnext11) && (tnext20 < tnext12) && (tnext20 < tnext13) &&(tnext20 < tnext14)  &&(tnext20 < tnext15) && (tnext20 < tnext16) &&(tnext20 < tnext17) && (tnext20 < tnext18) &&(tnext20 < tnext19)
        
        t = tnext20;        
        q20 = quantize(x20, q20, dQ20, eps20);
        x20 = dint(t, tlast20, x20, d20);
        d20 = f20(Rc, C,x10,x20);
        tnext20 = ta(d20, x20, q20, dQ20, t);
        
        
        tlast20 = t;
        if q20 ~= qlast20
            qlast20 = q20;
            thist20(k20) = t;
            qhist20(k20) = q20;
            k20 = k20 + 1;
            x10 = dint(t, tlast10, x10, d10);
            tlast10 =t;
            d10 = f10(Rs, L,x10, x19,x20);
            tnext10 = ta(d10, x10, q10, dQ10, t);
            t= min([tnext1,tnext2,tnext3,tnext4,tnext5,tnext6,tnext7,tnext8,tnext9,tnext10,tnext11,tnext12,tnext13,tnext14,tnext15,tnext16,tnext17,tnext18,tnext19,tnext20]);
     
        end
    elseif tnext1 == tnext2== tnext3== tnext4 == tnext5== tnext6== tnext7 == tnext8== tnext9== tnext10 == tnext11== tnext12== tnext13 == tnext14== tnext15== tnext16 == tnext17== tnext18== tnext19 == tnext20
        
        t = tnext1;
        x1 = dint(t, tlast1, x1, d1);
        d1 = f1(Rs, L, E, x1, x3);
        q1 = quantize(x1, q1, dQ1, eps1);
        tnext1 = ta(d1, x1, q1, dQ1, t);
        
        tlast1 = t;
        if q1 ~= qlast1
            qlast1 = q1;
            thist1(k1) = t;
            qhist1(k1) = q1;
            k1 = k1 + 1;
            tnext2 = t;
        end
        
        t = tnext2;
        x2 = dint(t, tlast2, x2, d2);
        d2 = f2(Rs, L,x2, x3,x4);
        q2 = quantize(x2, q2, dQ2, eps2);
        tnext2 = ta(d2, x2, q2, dQ2, t);
       
        tlast2 = t;
        if q2 ~= qlast2
            qlast2 = q2;
            thist2(k2) = t;
            qhist2(k2) = q2;
            k2 = k2 + 1;
            tnext1 = t;
        end
        
    end
    
end

% figure(3)
% plot(thist1, qhist1, 'm-');hold on;
plot(thist10, qhist10, 'g-');hold on;
% plot(thist3, qhist3, 'b-');hold on;
plot(thist20, qhist20, 'r-');
% legend('QSS_Current_x1','QSS_Current_x2','QSS_Voltage_x3','QSS_Voltage_x4')
legend('Matlab_Current_x10','Matlab_Voltage_x20','QSS_Current_x10','QSS_Voltage_x20')



function [x] = dint(t, tlast, x0, d)

    % update elapsed time for each component:
    dt = t - tlast;

    % update internal states:
    x = x0 + d * dt;
    
    
end

function [tnext] = ta(d, x, q, dQ, t)
    if (d > 0)
        dt=(q + dQ - x) / d;
        tnext = t + abs(dt);
                 
    elseif (d < 0)
        dt = (q - 0.5 * dQ - x) / d;
        tnext = t + abs(dt);
    else
        tnext = inf;
    end
    
     A=[d x q dQ t];
     
       if(tnext<0)
         disp(A);
%          display(d,x,q,dQ,t);
         error('tnext negavtive');
         
       end
      
       
end

%----------------------------------------------------------------------
function [d1] = f1(Rs, L, E, x1, x11)
    d1 = 1/L * (E - (x1*Rs) - x11);
end

function [d2] = f2(Rs, L,x2, x11,x12)
    d2 = 1/L * (-(x2*Rs) + x11-x12);
end

function [d3] = f3(Rs, L,x3, x12,x13)
    d3 = 1/L * (- (x3*Rs) + x12-x13);
end

function [d4] = f4(Rs, L,x4, x13,x14)
    d4 = 1/L * (- (x4*Rs) + x13-x14);
end

function [d5] = f5(Rs, L,x5, x14,x15)
    d5 = 1/L * (- (x5*Rs) + x14-x15);
end

function [d6] = f6(Rs, L,x6, x15,x16)
    d6 = 1/L * (- (x6*Rs) + x15-x16);
end

function [d7] = f7(Rs, L,x7, x16,x17)
    d7 = 1/L * (- (x7*Rs) + x16-x17);
end

function [d8] = f8(Rs, L,x8, x17,x18)
    d8 = 1/L * (- (x8*Rs) + x17-x18);
end

function [d9] = f9(Rs, L,x9, x18,x19)
    d9 = 1/L * (- (x9*Rs) + x18-x19);
end

function [d10] = f10(Rs, L,x10, x19,x20)
    d10 = 1/L * (- (x10*Rs) + x19-x20);
end
%----------------------------------------------------------------------
function [d11] = f11(Rc, C, x1,x2, x11)
    d11 = 1/C * (- ((x11)/Rc) + x1-x2);
end

function [d12] = f12(Rc, C,x2,x3,x12)
    d12 = 1/C * (- (x12/Rc) + x2-x3);
end

function [d13] = f13(Rc, C,x3,x4,x13)
    d13 = 1/C * (- (x13/Rc) + x3-x4);
end

function [d14] = f14(Rc, C,x4,x5,x14)
    d14 = 1/C * (- (x14/Rc) + x4-x5);
end

function [d15] = f15(Rc, C,x5,x6,x15)
    d15 = 1/C * (- (x15/Rc) + x5-x6);
end

function [d16] = f16(Rc, C,x6,x7,x16)
    d16 = 1/C * (- (x16/Rc) + x6-x7);
end

function [d17] = f17(Rc, C,x7,x8,x17)
    d17 = 1/C * (- (x17/Rc) + x7-x8);
end

function [d18] = f18(Rc, C,x8,x9,x18)
    d18 = 1/C * (- (x18/Rc) + x8-x9);
end

function [d19] = f19(Rc, C,x9,x10,x19)
    d19 = 1/C * (- (x19/Rc) + x9-x10);
end

function [d20] = f20(Rc, C,x10,x20)
    d20 = 1/C * (- (x20/Rc) + x10);
end

%-----------------------------------------------------------------------
function [q] = quantize(x, q, dQ, eps)

    if x >= q + dQ - eps
        q = q + dQ;
    elseif x <= q - 0.5 * dQ + eps
        q = q - dQ;
    end

end

