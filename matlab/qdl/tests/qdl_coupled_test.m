clear variables

dqmin = 0.1;
dqmax = 0.1;
dqerr = 0.001;

Va = 10;
Ra = 0.1;
La = 0.01;
Ke = 0.1;
Kt = 0.1;
Jm = 0.1;
Bm = 0.001;
Tm = 0.0;

sys = QdlSystem(dqmin, dqmax, dqerr);

va = QdlNode('Va (V) (100Hz d=0.5)', 0, 0, 0);
va.stim_type = QdlSystem.StimSource;
va.source_type = QdlSystem.SourcePWM;
va.freq = 100; 
va.duty = 0.5; 
va.v1 = Va; 
va.v2 = 0; 

gnd = QdlNode('gnd', 0, 0, 0);
gnd.stim_type = QdlSystem.StimSource;
gnd.source_type = QdlSystem.SourceDC;
gnd.vdc = 0.0; 

nmech = QdlNode('shaft speed (rad/s)', Jm, Bm, Tm);

barm = QdlBranch('armature current (A)', La, Ra, 0);
barm.connect(va, gnd);

barm.add_tnode(nmech, -Ke);
nmech.add_sbranch(barm, Kt);

sys.add_node(va);
sys.add_node(gnd);
sys.add_node(nmech);
sys.add_branch(barm);

sys.init();
sys.runto(10);
sys.H(3) = -0.5;
sys.runto(20);

nbins = 100;
ymax = 100;

figure;

subplot(3, 1, 1);
sys.plot(va, 0, 1, 1, 0, nbins, 0, 0, 0, ymax);

subplot(3, 1, 2);
sys.plot(barm, 0, 1, 1, 0, nbins, 0, 0, 0, ymax*10);

subplot(3, 1, 3);
sys.plot(nmech, 0, 1, 1, 0, nbins, 0, 0, 0, ymax);




