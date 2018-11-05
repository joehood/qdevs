

dQ = 0.01;
eps = dQ/4;

H = 1.0;
G = 1.0;
C = 1.0;
E = 1.0;
R = 1.0;
L = 1.0;

% create the system:
sys = QssLimSystem();

% define nodes:
gnd = QssLimNode(0, 0, 0, 0.0, 1, dQ, eps);
node1 = QssLimNode(H, G, C, 0.0, 0, dQ, eps);

% define branches:
branch1 = QssLimBranch(E, R, L, 0.0, node1, gnd, dQ, eps);

% add atoms to system:
sys.add_node(gnd);
sys.add_node(node1);
sys.add_branch(branch1);

% advance the simulator:
sys.init(0.0);

