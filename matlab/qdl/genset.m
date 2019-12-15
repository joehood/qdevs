clear variables

dqmin = 0.001;
dqmax = 0.001;
dqerr = 0.0001;

alpha = 0.5; 
alpha1 = 0.5; 
alpha2 = 2;
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

gnd = QdlNode('gnd', 0.00001, 0, 0);
%gnd.source_type = QdlSystem.SourceDC;
%gnd.vdc = 0.0;

GenNode = QdlNode('Generator Terminal Voltage (V)', C, G, 0);
MechNode = QdlNode('Generator Speed (rad/s)', C, G, 0);
GenBranch = QdlBranch('Generator Current (A)', Ls, Rs, Vs);
ConvBranch = QdlBranch('Converter 1 Primary Current (A)', L, R, 0);
ConvNode = QdlNode('Converter 1 Primary Voltage (V)', C, G, 0);
SecBranch1 = QdlBranch('Converter 1 Secondary Current (A)', L, R, 0);
PriNode2 = QdlNode('Converter 2 Primary Voltage (V)', C, G, 0);
SecBranch2 = QdlBranch('Converter 2 Secondary Current (A)', L, R, 0);
MotorNode = QdlNode('Motor Speed (rad/s)', J , B , 0);

% 3. controlled source coupling:

MechNode.add_sbranch(GenBranch, -alpha);
GenBranch.add_tnode(MechNode, alpha);
GenNode.add_sbranch(ConvBranch, -alpha);
ConvBranch.add_tnode(GenNode, alpha);
ConvNode.add_sbranch(SecBranch1, -alpha1);
SecBranch1.add_tnode(ConvNode, alpha1);
PriNode2.add_sbranch(SecBranch2, -alpha2);
SecBranch2.add_tnode(PriNode2, alpha2);

% 4. connect the system:

GenBranch.connect(gnd, GenNode);
ConvBranch.connect(gnd, ConvNode);
SecBranch1.connect(gnd, PriNode2);
SecBranch2.connect(gnd, MotorNode);

% 5. add atoms to system:

sys.add_node(gnd);
sys.add_node(MechNode);
sys.add_branch(GenBranch);
sys.add_node(GenNode);
sys.add_branch(ConvBranch);
sys.add_node(ConvNode);
sys.add_branch(SecBranch1);
sys.add_node(PriNode2);
sys.add_branch(SecBranch2);
sys.add_node(MotorNode);

sys.init();

% 6. qss simulation:

if 0
    sys.init();
    sys.runto(20);
    sys.E(1) = sys.E(1) * 1.2; % increase source voltage by 20%
    sys.runto(40);
end

% 7. state space simulation:

if 1
    dtss = 0.001;
    sys.E(1) = 10;
    sys.init_ss(dtss);
    [tss1, xss1] = sys.run_ss_to(dtss, 20);
    sys.E(1) = sys.E(1) * 1.2;
    [tss2, xss2] = sys.run_ss_to(dtss, 40);
end 

% 8. plot output:

figure;

% sys.plot(atom, dots, lines, upd, cumm_upd, bins, xlbl, tss, xss, ymax)

upd = 1;
cumm_upd = 0;
bins = 500;

subplot(4, 2, 1);
sys.plot(MechNode, 0, 1, upd, cumm_upd, bins, 0, [tss1(1,:), tss2(1,:)], [xss1(1,:), xss2(1,:)], 0);

subplot(4, 2, 2);
sys.plot(GenBranch, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

subplot(4, 2, 3);
sys.plot(GenNode, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

subplot(4, 2, 4);
sys.plot(ConvBranch, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

subplot(4, 2, 5);
sys.plot(ConvNode, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

subplot(4, 2, 6);
sys.plot(SecBranch1, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

subplot(4, 2, 7);
sys.plot(SecBranch2, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

subplot(4, 2, 8);
sys.plot(MotorNode, 0, 1, upd, cumm_upd, bins, 0, 0, 0, 0);

figure;

plot(tss1(1,:), xss1(1,:))



