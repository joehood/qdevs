clear variables

dqmin = 0.001;
dqmax = 0.001;
dqerr = 0.001;

sys = QdlSystem(dqmin, dqmax, dqerr);

nn = 11; % includes ground
nb = 10;
R = 1;
E1 = 5; E2 = 10;
G = 1/10;
C = [0.09, 0.09, 0.09 ,0.09, 0.09, 0.09, 0.09, 0.01, 0.01, 0.9];
L = [0.001, 0.001, 0.001, 0.0001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01];

nodes = QdlNode.empty(0);
branches = QdlBranch.empty(0);

nodes(1) = QdlNode('gnd', 100000, 100000, 0);
sys.add_node(nodes(1));

for k = 2:nn
    nodes(k) = QdlNode(strcat('node', num2str(k-1)), C(k-1), G, 0);
    sys.add_node(nodes(k));
end

for k = 1:nb
    branches(k) = QdlBranch(strcat('node', num2str(k)), L(k), R, 0);
    sys.add_branch(branches(k));
    branches(k).connect(nodes(k), nodes(k+1));
end

sys.init();
sys.E(1) = E1;
[t1, x1] = sys.run_ss(0.01, 20);
sys.runto(20);
sys.E(1) = E2;
[t2, x2] = sys.run_ss(0.01, 40);
sys.runto(40);

figure;

sys.plot(nodes(5), 1, 0, 1, 0, 500, 0, [t1,t2], [x1(5,:),x2(5,:)], 0);

%for k = 1:nn
%    sys.plot(nodes(5), 1, 0, 1, 0, 500, 0, 0, 0, 0);
%end






