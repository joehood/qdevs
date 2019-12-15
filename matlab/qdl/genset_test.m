clear variables

dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.001;

Rs = 0.01;
Ld = 0.001;
Lq = 0.001;
Rf = 0.01;
Lf = 0.001;
Cf = 0.001;
J = 10;
P = 2;
B = 0.01;
Ed = 100;
Eq = 0;
Rl = 1;
w0 = 60;
Pm = 1e3;
kt = 0.5;
ke = 0.5;

sys = QdlSystem(dqmin, dqmax, dqerr);

gnd = QdlNode('gnd', 0, 0, 0);
gnd.source_type = QdlSystem.SourceDC;
gnd.vdc = 0.0;

Dbranch = QdlBranch('Dbranch', Ld, Rs, Ed);
ConvBranch = QdlBranch('Convbranch', Lf, Rf, 0);
GenNode = QdlNode('GenNode', 0.001, 0, 0);
SinkNode = QdlNode('SinkNode', Cf, Rl, 0);
MechNode = QdlNode('MechNode', 2*J/P, B, Pm/w0);

% 3. controlled sources:

MechNode.add_sbranch(Dbranch, -kt);
Dbranch.add_tnode(MechNode, -ke);

% 4. wiring of the system:

Dbranch.connect(gnd, GenNode);
ConvBranch.connect(gnd, SinkNode);

% 5. add atoms to system:

sys.add_node(gnd);
sys.add_branch(Dbranch);
sys.add_branch(ConvBranch);
sys.add_node(GenNode);
sys.add_node(SinkNode);
sys.add_node(MechNode);

sys.init();
sys.runto(0.01);

figure;


subplot(4, 2, 1);
sys.plot(Dbranch, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 2);
sys.plot(GenNode, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 3);
sys.plot(ConvBranch, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 4);
sys.plot(SinkNode, 0, 1, 0, 0, 0, 0, 0, 0, 0);




