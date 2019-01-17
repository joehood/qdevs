 clear variables

dqmin = 0.001;
dqmax = 0.001;
dqerr = 0.001;

sys = QdlSystem(dqmin, dqmax, dqerr);

n1 = QdlNode('Node 1 Voltage (V)', 1, 1, 1);
n2 = QdlNode('Node 2 Voltage (V)', 1, 1, 1);
n3 = QdlNode('Node 3 Voltage (V)', 1, 1, 1);
n4 = QdlNode('Node 4 Voltage (V)', 1, 1, 1);
n5 = QdlNode('Node 5 Voltage (V)', 1, 1, 1);
b1 = QdlBranch('Branch 1-2 Current (A)', 1, 2, 1);
b2 = QdlBranch('Branch 2-3 Current (A)', 1, 1.5, 1);
b3 = QdlBranch('Branch 3-4 Current (A)', 1, 1, 1);
b4 = QdlBranch('Branch 4-5 Current (A)', 1, .5, 1);

b1.connect(n1, n2);
b2.connect(n2, n3);
b3.connect(n3, n4);
b4.connect(n4, n5);

sys.add_node(n1);
sys.add_node(n2);
sys.add_node(n3);
sys.add_node(n4);
sys.add_node(n5);
sys.add_branch(b1);
sys.add_branch(b2);
sys.add_branch(b3);
sys.add_branch(b4);

sys.init();
[t1, x1] = sys.run_ss(0.01, 10);
sys.runto(10);
sys.E(1) = 2;
[t2, x2] = sys.run_ss(0.01, 200);
sys.runto(20);


figure;
r = 5;
c = 1;

subplot(r, c, 1);
sys.plot(n1, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(1,:),x2(1,:)], 0);

subplot(r, c, 2);
sys.plot(n2, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(2,:), x2(2,:)], 0);

subplot(r, c, 3);
sys.plot(n3, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(3,:), x2(3,:)], 0);

subplot(r, c, 4);
sys.plot(n4, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(4,:), x2(4,:)], 0);

subplot(r, c, 5);
sys.plot(n5, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(5,:), x2(5,:)], 0);

figure;
r = 4;
c = 1;

subplot(r, c, 1);
sys.plot(b1, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(6,:), x2(6,:)], 0);

subplot(r, c, 2);
sys.plot(b2, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(7,:), x2(7,:)], 0);

subplot(r, c, 3);
sys.plot(b3, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(8,:), x2(8,:)], 0);

subplot(r, c, 4);
sys.plot(b4, 1, 0, 1, 0, 400, 0, [t1, t2], [x1(9,:), x2(9,:)], 0);





