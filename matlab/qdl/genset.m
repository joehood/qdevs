clear variables

dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.001;

alpha = 0.5; 
alpha1 = 0.5; 
alpha2=2;

R = 0.1;
C = 0.001;
L = 0.01;
G = 0.1;

B = 1.0;
J = 10.0; 

Vs = 10;
Rs = 0.001;
Ls = 0.01;

% 1. define the system:
sys = QdlSystem(dqmin, dqmax, dqerr);

% 2. define the atoms:

gnd = QdlNode('gnd', 0, 0, 0);
gnd.source_type = QdlSystem.SourceDC;
gnd.vdc = 0.0;
%----------------------------------------------------------------------
MechNode = QdlNode('Generator Input Voltage (V)', C, G, 0);


GenBranch = QdlBranch('Source Current (A)', Ls, Rs, Vs);
GenNode = QdlNode('Generator Primary Voltage (V)', C, G, 0);
ConvBranch = QdlBranch('Generator Secondary Current (A)', L, R, 0);

ConvNode = QdlNode('Generator Output Voltage (V)', C, G, 0);

%first converter

%PriNode1 = QdlNode('Converter 1 Primary Voltage (V)', C, G, 0);
SecBranch1 = QdlBranch('Converter 1 Secondary Current (A)', L, R, 0);

%------------------------------------------------------------------------
%Second converter

PriNode2 = QdlNode('Converter 2 Primary Voltage (V)', C, G, 0);
SecBranch2 = QdlBranch('Converter 2 Secondary Current (A)', L, R, 0);

%-------------------------------------------------------------------------

MotorNode = QdlNode('Motor Speed (rad/s)',J , B , 0);


% 3. controlled sources:

MechNode.add_sbranch(GenBranch, -alpha);
GenBranch.add_tnode(MechNode, alpha);

GenNode.add_sbranch(ConvBranch, -alpha);
ConvBranch.add_tnode(GenNode, alpha);


ConvNode.add_sbranch(SecBranch1, -alpha1);
SecBranch1.add_tnode(ConvNode, alpha1);
%-------------------------------
PriNode2.add_sbranch(SecBranch2, -alpha2);
SecBranch2.add_tnode(PriNode2, alpha2);

% 4. wiring of the system:

GenBranch.connect(gnd, GenNode);
ConvBranch.connect(gnd, ConvNode);
%SourceBranch.connect(gnd, PriNode1);
SecBranch1.connect(gnd, PriNode2);
SecBranch2.connect(gnd,MotorNode);



% 5. add atoms to system:


sys.add_node(gnd);
sys.add_node(MechNode);
sys.add_branch(GenBranch);
sys.add_node(GenNode);
sys.add_branch(ConvBranch);
sys.add_node(ConvNode);
sys.add_node(gnd);

%sys.add_branch(SourceBranch);

sys.add_node(ConvNode);
sys.add_branch(SecBranch1);

sys.add_node(PriNode2);
sys.add_branch(SecBranch2);

sys.add_node(MotorNode);

sys.init();
sys.runto(10);
sys.E(1) = 10 * 1.1; % increase source voltage
sys.runto(20);

figure;

subplot(4, 2, 1);
sys.plot(GenBranch, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 2);
sys.plot(GenNode, 0, 1, 0, 0, 0, 0, 0, 0, 0);


subplot(4, 2, 3);
sys.plot(ConvBranch, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 4);
sys.plot(ConvNode, 0, 1, 0, 0, 0, 0, 0, 0, 0);

%subplot(4, 2, 5);
%sys.plot(PriNode1, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 5);
sys.plot(SecBranch1, 0, 1, 0, 0, 0, 0, 0, 0, 0);


subplot(4, 2, 6);
sys.plot(PriNode2, 0, 1, 0, 0, 0, 0, 0, 0, 0);

subplot(4, 2, 7);
sys.plot(SecBranch2, 0, 1, 0, 0, 0, 0, 0, 0, 0);





