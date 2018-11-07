clear all

dQ = 0.005;

tstop = 15.0;

H = 1.0;
G = 0.1;
C = 0.01;
E = 5.0;
R = 1.0;
L = 0.01;

h = 0.001;

n = 10;

a11 = eye(n) * -R/L;
a12 = eye(n) * -1/L;
a21 = eye(n) *  1/C;
a22 = eye(n) * -G/C;

for i = 1:n-1
   a21(i, i+1) = -1/C;
end

for i = 2:n
   a12(i, i-1) =  1/L;
end

a = [a11, a12 ; a21, a22];
 
b = [ 1/L, zeros(1, 2*n-1)]';

u = [ E ]';

apr = inv(eye(20) - a*h);
bpr = apr* h * b ;


t = 0:h:tstop;
x = zeros(20, length(t));

for k = 2:length(t)
    if t(k) > 8
        u = 10;
    end
    x(:,k) = apr * x(:,k-1) + bpr * u;
    
end

% create the system:
sys = QssLimSystem(dQ, dQ, 0.25);

% define nodes:
gnd = QssLimGround(0);
nodes = [ QssLimNode('node1', 0, G, C, 0), 
          QssLimNode('node2', 0, G, C, 0),
          QssLimNode('node3', 0, G, C, 0),
          QssLimNode('node4', 0, G, C, 0),
          QssLimNode('node5', 0, G, C, 0),
          QssLimNode('node6', 0, G, C, 0),
          QssLimNode('node7', 0, G, C, 0),
          QssLimNode('node8', 0, G, C, 0),
          QssLimNode('node9', 0, G, C, 0),
          QssLimNode('node10', 0, G, C, 0) ];

% define and connect branches:
branches = [ QssLimBranch('branch1', E, R, L, 0, gnd, nodes(1)), 
             QssLimBranch('branch2', 0, R, L, 0, nodes(1), nodes(2)), 
             QssLimBranch('branch3', 0, R, L, 0, nodes(2), nodes(3)), 
             QssLimBranch('branch4', 0, R, L, 0, nodes(3), nodes(4)), 
             QssLimBranch('branch5', 0, R, L, 0, nodes(4), nodes(5)), 
             QssLimBranch('branch6', 0, R, L, 0, nodes(5), nodes(6)), 
             QssLimBranch('branch7', 0, R, L, 0, nodes(6), nodes(7)), 
             QssLimBranch('branch8', 0, R, L, 0, nodes(7), nodes(8)), 
             QssLimBranch('branch9', 0, R, L, 0, nodes(8), nodes(9)), 
             QssLimBranch('branch10', 0, R, L, 0, nodes(9), nodes(10)) ];

% add atoms to system:
sys.add_ground(gnd);
for k = 1:10
    sys.add_node(nodes(k));
    sys.add_branch(branches(k));
end

% initialize and run:
sys.init(0);
sys.run(8.0);
branches(1).E = 10.0;
sys.run(tstop);

figure;

plot(t, x(20,:), 'b-'); hold on;
plot(nodes(10).thist, nodes(10).qhist, 'k.');

disp('done.');

function plot_qss(thist, qhist, t, linestyle, pointstyle)
    v1 = interp1(thist, qhist, t, 'previous');
    plot(t, v1, linestyle); hold on;
    plot(thist, qhist, pointstyle);
end


