clear all

R = 1;
L = 1;
E = 1;
G = 1;
C = 1;
H = 1;

t = 0.0;
tstop = 10;

dQ = 0.01;
eps = dQ/4;

%node = QssOdeAtom(dQ, eps, C, G, H);
%branch = QssOdeAtom(dQ, eps, L, R, E);
node = LiqssOdeAtom(dQ, eps, C, G, H);
branch = LiqssOdeAtom(dQ, eps, L, R, E);

node.init();
branch.init();

node.update(0, tstop);
branch.update(0, tstop);

while (t < tstop)

    % determine and cache inputs for atoms from system connectivity:
    vnode = node.q;
    ibranch = -branch.q;
    
    tnext_node = node.tnext;
    tnext_branch = branch.tnext;
    
    tnext = min([tnext_node, tnext_branch]);
    t = tnext;
 
    if tnext_node == tnext
        if node.update(ibranch, tstop) branch.tnext = t; end  
    end
    
    if tnext_branch == tnext      
        if branch.update(vnode, tstop) node.tnext = t; end 
    end
    
end

% add final value of history at tstop time:
node.thist(node.k) = tstop;
node.qhist(node.k) = node.qhist(node.k-1);
branch.thist(branch.k) = tstop;
branch.qhist(branch.k) = branch.qhist(branch.k-1);

% state space reference simulation:
h = 0.01;
A = [-R/L  -1/L ; 1/C  -G/C ];
b = [ 1/L 0 ; 0 1/C];
u = [ E  H ]';
Apr = inv(eye(2) - h * A);
bpr = h * b * Apr;
t = 0:h:tstop;
n = length(t);
x = zeros(2, n);  % [ vnode  ibranch ]
for k=2:n
    x(:,k) = Apr*x(:,k-1) + bpr*u; 
end

t2 = 0:0.001:tstop;
 
plot_qss(node.thist, node.qhist, t2, 'r-', 'k.');
plot_qss(branch.thist, branch.qhist, t2, 'b-', 'k.');

plot(t, x(1,:), 'm--'); 
plot(t, x(2,:), 'c--'); 

disp('done.')

function plot_qss(thist, qhist, t, linestyle, pointstyle)
    v1 = interp1(thist, qhist, t, 'previous');
    plot(thist, qhist, pointstyle); hold on;
    plot(t, v1, linestyle);
end
