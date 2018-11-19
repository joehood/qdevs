clear all

dQVoltage = 0.0001;
dQCurrent = 0.0001;
eps_factor = 0.25;
epsVoltage = eps_factor * dQVoltage;
epsCurrent = eps_factor * dQCurrent;

tstop = 100.0;

H = 1.03;
G = 0.1;
C = 0.013;
E = 5.0;
R = 1.0;
L = 0.012;

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
sys = LiqssLimSystem('20 order system');

% define nodes:

gnd = LiqssLimGround('grd', 0);

nodes = [ LiqssLimNode('node1',  0, G, C, 0, dQVoltage, epsVoltage), 
          LiqssLimNode('node2',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node3',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node4',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node5',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node6',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node7',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node8',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node9',  0, G, C, 0, dQVoltage, epsVoltage),
          LiqssLimNode('node10', 0, G, C, 0, dQVoltage, epsVoltage) ];

% define and connect branches:

branches = [ LiqssLimBranch('branch1',  E, R, L, 0, gnd,      nodes(1),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch2',  0, R, L, 0, nodes(1), nodes(2),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch3',  0, R, L, 0, nodes(2), nodes(3),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch4',  0, R, L, 0, nodes(3), nodes(4),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch5',  0, R, L, 0, nodes(4), nodes(5),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch6',  0, R, L, 0, nodes(5), nodes(6),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch7',  0, R, L, 0, nodes(6), nodes(7),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch8',  0, R, L, 0, nodes(7), nodes(8),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch9',  0, R, L, 0, nodes(8), nodes(9),  dQCurrent, epsCurrent), 
             LiqssLimBranch('branch10', 0, R, L, 0, nodes(9), nodes(10), dQCurrent, epsCurrent) ];

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

rows = 2
cols = 2

subplot(rows, cols, 1);
plot_qss(branches(1).thist, branches(1).qhist, x(1, :), t, 'branch (A)', 'Branch 1');

subplot(rows, cols, 2);
plot_qss(branches(10).thist, branches(10).qhist, x(10, :), t, 'branch (A)', 'Branch 10');

subplot(rows, cols, 3);
plot_qss(nodes(1).thist, nodes(1).qhist, x(11, :), t, 'voltage (V)', 'Node 1');

subplot(rows, cols, 4);
plot_qss(nodes(10).thist, nodes(10).qhist, x(20, :), t, 'voltage (V)', 'Node 10');

disp('done.');

function plot_qss(thist, qhist, xss, t, label, name)

    upd_max = 100
    binwidth = 0.01

    yyaxis left
    
    plot(t, xss, 'c--'); 
    
    hold on
    
    v1 = interp1(thist, qhist, t, 'previous');
    plot(thist, qhist, 'k.');
    plot(t, v1, 'b-');
    
    ylabel(label);

    yyaxis right
    
    h1 = histogram(thist);
    h1.BinWidth = binwidth;
    h1.LineStyle = 'none';
    
    ax = gca;
    %ax.XLim = [0 1000];
    %ax.YLim = [0 upd_max];
    
    ylabel(strcat('updates per', {'  '}, num2str(binwidth), ' s'));

    xlabel('t (s)');
    
    legend({'ss', 'qss', 'qss (zoh)', 'updates'}, 'Location', 'east');
    
    title(name);

end


