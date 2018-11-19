clear all

%
%                               +
%          .-------UUU--VVV--( )------.
%          |                          |
%     .----+----.                .----+----.
%     |    |    |                |    |    |
%    (v)  [ ]  ===              (v)  [ ]  ===
%     |    |    |                |    |    |
%     '----+----'                '----+----'
%         ---                        ---
%          '                          '
%

R = 1;
L = 1;
E = 1;

G = 1;
C1 = 1;
C2 = 100;
H = 1;

t = 0.0;
tstop = 52;

dQ = 0.01;
eps = dQ/4;

if 0
    node1 = QssOdeAtom(dQ, eps, C1, G, H);
    node2 = QssOdeAtom(dQ, eps, C2, G, H);
    branch12 = QssOdeAtom(dQ, eps, L, R, E);
else
    node1 = LiqssOdeAtom(dQ, eps, C1, G, H);
    node2 = LiqssOdeAtom(dQ, eps, C2, G, H);
    branch12 = LiqssOdeAtom(dQ, eps, L, R, E);
end

node1.init();
node2.init();
branch12.init();

node1.update(0, tstop);
node2.update(0, tstop);
branch12.update(0, tstop);

while (t < tstop)

    % determine and cache inputs for atoms from system connectivity:
    v12 = node1.q - node2.q;
    isum1 = -branch12.q;
    isum2 = branch12.q;
    
    tnext1 = node1.tnext;
    tnext2 = node2.tnext;
    tnext12 = branch12.tnext;
    
    tnext = min([tnext1, tnext2, tnext12]);
    
    t = tnext
 
    if tnext1 == tnext     
        if node1.update(isum1, tstop) branch12.tnext = t; end  
    end
    
    if tnext2 == tnext       
        if node2.update(isum2, tstop) branch12.tnext = t; end  
    end
    
    if tnext12 == tnext       
        if branch12.update(v12, tstop)
            node1.tnext = t;
            node2.tnext = t;
        end 
    end
    
end

% add final value of history at tstop time:
node1.thist(node1.k) = tstop;
node1.qhist(node1.k) = node1.qhist(node1.k-1);
node2.thist(node2.k) = tstop;
node2.qhist(node2.k) = node2.qhist(node2.k-1);
branch12.thist(branch12.k) = tstop;
branch12.qhist(branch12.k) = branch12.qhist(branch12.k-1);

% state space reference simulation:
h = 0.01;

A = [ -G/C1   1/C1   0   ;
      -1/L  -R/L   1/L ;
       0    -1/C2  -G/C2 ];                 
b = [  1/C1   0     0   ;
       0     1/L   0   ;
       0     0     1/C2 ];
   
u = [ E  H  E ]';

Apr = inv(eye(3) - h * A);
bpr = h * b * Apr;

t = 0:h:tstop;
n = length(t);
x = zeros(3, n);  % [ vnode1  ibranch12  vnode2 ]

for k=2:n
    x(:,k) = Apr*x(:,k-1) + bpr*u; 
end


t2 = 0:0.001:tstop;

plot_qss(node1.thist, node1.qhist, t2, 'r-', 'k.');
plot_qss(node2.thist, node2.qhist, t2, 'b-', 'k.');
plot_qss(branch12.thist, branch12.qhist, t2, 'g-', 'k.');

plot(t, x(1,:), 'c--');
plot(t, x(2,:), 'k--');
plot(t, x(3,:), 'm--');

disp('done.');

function plot_qss(thist, qhist, t, linestyle, pointstyle)
    %v1 = interp1(thist, qhist, t, 'previous');
    plot(thist, qhist, pointstyle); hold on;
    %plot(t, v1, linestyle);
end

