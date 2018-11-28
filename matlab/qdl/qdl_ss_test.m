clear variables

dqmin = 0.001;
dqmax = 0.001;
dqerr = 0.001;

sys = QdlSystem(dqmin, dqmax, dqerr);

n1 = QdlNode('n1', 1, 1, 1);
n2 = QdlNode('n2', 1, 1, 1);
n3 = QdlNode('n3', 1, 1, 1);
n4 = QdlNode('n4', 1, 1, 1);
n5 = QdlNode('n5', 1, 1, 1);
b1 = QdlBranch('b1', 1, 2, 1);
b1.connect(n1, n2);
b2 = QdlBranch('b2', 1, 1.5, 1);
b2.connect(n2, n3);
b3 = QdlBranch('b3', 1, 1, 1);
b3.connect(n3, n4);
b4 = QdlBranch('b4', 1, .5, 1);
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
sys.runto(10);
%sys.E(1) = 2;
%sys.runto(200);

[t, x] = sys.run_ss(0.01, 10);

figure;

r = 5;
c = 2;

subplot(r, c, 1);
sys.plot(n1, 1, 1, 0, 500, 0, t, x(1,:), 0);

subplot(r, c, 3);
sys.plot(n2, 1, 1, 0, 500, 0, t, x(2,:), 0);

subplot(r, c, 5);
sys.plot(n3, 1, 1, 0, 500, 0, t, x(3,:), 0);

subplot(r, c, 7);
sys.plot(n4, 1, 1, 0, 500, 0, t, x(4,:), 0);

subplot(r, c, 9);
sys.plot(n5, 1, 1, 0, 500, 0, t, x(5,:), 0);

subplot(r, c, 2);
sys.plot(b1, 1, 1, 0, 500, 0, t, x(6,:), 0);

subplot(r, c, 4);
sys.plot(b2, 1, 1, 0, 500, 0, t, x(7,:), 0);

subplot(r, c, 6);
sys.plot(b3, 1, 1, 0, 500, 0, t, x(8,:), 0);

subplot(r, c, 8);
sys.plot(b4, 1, 1, 0, 500, 0, t, x(9,:), 0);





