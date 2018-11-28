clear variables

dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.001;

Vs = 10;
Gs = 1e2;
Clim = 1e-6;
Ra = 0.1;
La = 0.001;
Ke = 0.1;
Kt = 0.1;
Jm = 0.01;
Bm = 0.001;
Jp = 0.01;
Fp = 1;
JL = 0.5
TL = 0;
BL = 0.1;

sys = QdlSystem(dqmin, dqmax, dqerr);

vs = QdlNode('terminal', Clim, Gs, Vs*Gs);
vs.source_type = QdlSystem.SourceDC;
vs.vdc = 10.0;

vm = QdlNode('rotor', Jm, Bm, 0);
vl = QdlNode('load', JL, BL, TL);

vg = QdlNode('ground', 0, 0, 0);
vg.source_type = QdlSystem.SourceDC;
vg.vdc = 0;

ba = QdlBranch('armature', La, Ra, 0);
bs = QdlBranch('shaft', Jp, Fp, 0);

ba.connect(vs, vg);
bs.connect(vm, vl);

ba.add_tnode(vm, -Ke);
vm.add_sbranch(ba, Kt);

sys.add_node(vs);
sys.add_node(vm);
sys.add_node(vl);
sys.add_node(vg);
sys.add_branch(ba);
sys.add_branch(bs);

sys.init();
sys.runto(20);
sys.H(2) = -2;
sys.runto(30);

figure;

r = 3;
c = 2;
nbins = 1000;
ymax = 100;

subplot(r, c, 1);
sys.plot(vs, 0, 1, 0, nbins, 0, 0, 0, ymax);

subplot(r, c, 2);
sys.plot(ba, 0, 1, 0, nbins, 0, 0, 0, ymax);

subplot(r, c, 3);
sys.plot(vm, 0, 1, 0, nbins, 0, 0, 0, ymax);

subplot(r, c, 4);
sys.plot(bs, 0, 1, 0, nbins, 0, 0, 0, ymax);

subplot(r, c, 5);
sys.plot(vl, 0, 1, 0, nbins, 0, 0, 0, ymax);



