clear variables

dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.001;

alpha = 0.5; 

R = 0.1;
C = 0.001;
L = 0.01;
G = 0.1;

B = 1.0;
J = 10.0; 

Vs = 10;
Rs = 0.001;
Ls = 0.01;

sys = QdlSystem(dqmin, dqmax, dqerr);

gnd = QdlNode('gnd', 0, 0, 0);
gnd.source_type = QdlSystem.SourceDC;
gnd.vdc = 0.0; 
 
PriNode = QdlNode('Converter Primary Voltage (V)', C, G, 0);
SecBranch = QdlBranch('Converter Secondary Current (A)', L, R, 0);
SourceBranch = QdlBranch('Converter Source Current (A)', Ls, Rs, Vs);
MotorNode = QdlNode('Motor Speed (rad/s)',J , B , 0);

SecBranch.add_tnode(PriNode, alpha);
PriNode.add_sbranch(SecBranch, -alpha);

SourceBranch.connect(gnd,PriNode);
SecBranch.connect(gnd,MotorNode);

sys.add_node(PriNode);
sys.add_node(MotorNode);
sys.add_node(gnd);
sys.add_branch(SecBranch);
sys.add_branch(SourceBranch);

sys.init();
sys.runto(10);

figure;

subplot(2, 2, 1);
sys.plot(PriNode, 0, 1, 1, 0, 0, 0, 0, 0, 0);

subplot(2, 2, 2);
sys.plot(MotorNode, 0, 1, 1, 0, 0, 0, 0, 0, 0);


subplot(2, 2, 3);
sys.plot(SourceBranch, 0, 1, 1, 0, 0, 0, 0, 0, 0);

subplot(2, 2, 4);
sys.plot(SecBranch, 0, 1, 1, 0, 0, 0, 0, 0, 0);







