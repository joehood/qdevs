clear all


%
%          R     L     v1    R     L     v2
%    .----VVV---UUU----o----VVV---UUU----o
%    |      --->       |       --->      |
%   ( ) E    i1       [ ] G,C   i2      [ ] G,C
%    |                 |                 |
%    '-----------------+-----------------'
%                     -+- 
%                      '
%

dQvoltage = 0.01;
dQcurrent = 0.008;
eps_factor = 0.25;

tstop = 15.0;

H = 0.0;
G = 1.0;
C = 1.0;
E = 1.0;
R = 1.0;
L = 1.0;

% create the system:
sys = QssLimSystem(dQvoltage, dQcurrent, eps_factor);

% define nodes:
gnd = QssLimGround(0);
node1 = QssLimNode('node1', 0, G, C, 0);
node2 = QssLimNode('node2', 0, G, C, 0);

% define and connect branches:
branch1 = QssLimBranch('branch1', E, R, L, 0, gnd, node1);
branch2 = QssLimBranch('branch2', 0, R, L, 0, node1, node2);

% add atoms to system:
sys.add_ground(gnd);
sys.add_node(node1);
sys.add_node(node2);
sys.add_branch(branch1);
sys.add_branch(branch2);

% initialize and run:
sys.init(0);
sys.run(tstop*0.5);
branch1.E = 2.0;
sys.run(tstop);

%   state space benchmark:
%
%   E  = i1*R + i1'*L + v1
%   v1 = i2*R + i2'*L + v2
%   i1 = v1*G + v1'*C + i2
%   i2 = v2*G + v2'*C
%
%   i1' = 1/L * (E - i1*R - v1)
%   i2' = 1/L * (v1 - i2*R - v2)
%   v1' = 1/C * (i1 - v1*G - i2)
%   v2' = 1/C * (i2 - v2*G)
%
%   i1' = i1 * -R/L  +  i2 *  0    +  v1 * -1/L  +  v2 *  0
%   i2' = i1 *  0    +  i2 * -R/L  +  v1 *  1/L  +  v2 * -1/L
%   v1' = i1 *  1/C  +  i2 * -1/C  +  v1 * -G/C  +  v2 *  0
%   v2' = i1 *  0    +  i2 *  1/C  +  v1 *  0    +  v2 * -G/C

n = 4;
h = 0.001;
t = 0:h:tstop;

a = [ -R/L   0    -1/L   0   ;
       0    -R/L   1/L  -1/L ;
       1/C  -1/C  -G/C   0   ;
       0     1/C   0    -G/C ];
   
b = [  1/L   0     0     0   ]';

u = [  E  ]';

apr = inv(eye(n) - a * h);
bpr = apr * h * b;

x = zeros(n, length(t));

for k = 2:length(t)
    if t(k-1) > tstop * 0.5;
        u = [ 2.0 ];
    end
    x(:,k) = apr * x(:,k-1) + bpr * u; 
end

figure;

subplot(2, 2, 1);
plot(t, x(3,:), 'm--'); hold on;
plot_qss(node1.thist, node1.qhist, t, 'r-', 'k.');
legend({'ss', 'qss', 'qss (zoh)'}, 'Location', 'southeast');
title(node1.name);

subplot(2, 2, 2);
plot(t, x(4,:), 'm--'); hold on;
plot_qss(node2.thist, node2.qhist, t, 'r-', 'k.');
legend({'ss', 'qss', 'qss (zoh)'}, 'Location', 'southeast');
title(node2.name);

subplot(2, 2, 3);
plot(t, x(1,:), 'c--'); hold on;
plot_qss(branch1.thist, branch1.qhist, t, 'b-', 'k.');
legend({'ss', 'qss', 'qss (zoh)'}, 'Location', 'southeast');
title(branch1.name);

subplot(2, 2, 4);
plot(t, x(2,:), 'c--'); hold on;
plot_qss(branch2.thist, branch2.qhist, t, 'b-', 'k.');
legend({'ss', 'qss', 'qss (zoh)'}, 'Location', 'southeast');
title(branch2.name);


disp('done.');


function plot_qss(thist, qhist, t, linestyle, pointstyle)
    v1 = interp1(thist, qhist, t, 'previous');
    plot(thist, qhist, pointstyle);
    plot(t, v1, linestyle); hold on;
    
end




