clear variables

dqmin = .01;
dqmax = .01;
dqerr = 0.001;

nnodes = 5; % must be atleast 2
nodes = QdlNode.empty(0);
branches = QdlBranch.empty(0);

sys = QdlSystem(dqmin, dqmax, dqerr);

nodes(1) = QdlNode('node0', 0, 0, 0);

% dc:
if 1
nodes(1).source_type = QdlSystem.SourceDC;
nodes(1).vdc = 10.0;
end   

% pwm:
if 0
nodes(1).source_type = QdlSystem.SourcePWM;
nodes(1).freq = 10000;
nodes(1).duty = 0.5;
nodes(1).v1 = 0; 
nodes(1).v2 = 1;
end

% sine:
if 0
nodes(1).source_type = QdlSystem.SourceSINE;
nodes(1).freq = 10;
nodes(1).va = 1;
end

sys.add_node(nodes(1));

C = 1e-6;
L = 1e-6;
s = 2;

for k=2:nnodes-1
    nn = strcat('node', num2str(k-1));
    nodes(k) = QdlNode(nn, C*10^(k*s), 1, 0);
    nb = strcat('branch', num2str(k-1), '-', num2str(k));
    branches(k-1) = QdlBranch(nb, L*10^(k*s), 0.01, 0);
    branches(k-1).connect(nodes(k-1), nodes(k));
    sys.add_node(nodes(k));
    sys.add_branch(branches(k-1));
end

nodes(nnodes) = QdlNode(strcat('node', num2str(nnodes-1)), C*10^(nnodes*s), 1, 0);
nb = strcat('branch', num2str(nnodes-1), '-', num2str(nnodes));
branches(nnodes-1) = QdlBranch(nb, L*10^(nnodes*s), 0.01, 0);
branches(nnodes-1).connect(nodes(nnodes-1), nodes(nnodes));
sys.add_node(nodes(nnodes));
sys.add_branch(branches(nnodes-1));  

C*10^(nnodes*s)
L*10^(nnodes*s)

sys.init();
sys.runto(5000);
sys.xdc(1) = 20;
sys.runto(10000);

figure;


subplot(4, 2, 1);
sys.plot(nodes(2), 0, 1, 1, 0, 1000, 0, 0, 0, 0);
xlabel('t (s)');

subplot(4, 2, 2);
sys.plot(nodes(2), 0, 1, 1, 0, 10000000, 0, 0, 0, 0);
xlabel('t (s)');

subplot(4, 2, 3);
sys.plot(nodes(3), 0, 1, 1, 0, 1000, 0, 0, 0, 0);
xlabel('t (s)');

subplot(4, 2, 4);
sys.plot(nodes(3), 0, 1, 1, 0, 100000, 0, 0, 0, 0);
xlabel('t (s)');

subplot(4, 2, 5);
sys.plot(nodes(4), 0, 1, 1, 0, 1000, 0, 0, 0, 0);
xlabel('t (s)');

subplot(4, 2, 6);
sys.plot(nodes(4), 0, 1, 1, 0, 1000, 0, 0, 0, 0);
xlabel('t (s)');

subplot(4, 2, 7);
sys.plot(nodes(5), 0, 1, 1, 0, 100, 0, 0, 0, 0);
xlabel('t (s)');





