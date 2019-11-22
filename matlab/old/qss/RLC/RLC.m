clear variables

dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.001;

Vs = 10;
R = 5;
L = 1e-3;
C = 1e-3;

sys=QdlSystem(dqmin,dqmax,dqerr);

gnd = QdlNode('ground', 0, 0, 0);
gnd.source_type = QdlSystem.SourceDC;
gnd.vdc = 0;


vc=QdlNode('Capacitor node',C,0,0);
vc.source_type=QdlSystem.SourceDC;
vc.vdc=0;

br=QdlBranch('branch RL',L,R,Vs);



br.connect(gnd,vc);

sys.add_node(vc);
sys.add_branch(br);


sys.init();
sys.runto(10);


figure;

r = 3;
c = 2;
nbins = 1000;
ymax = 100;
dots = 0;

sys.plot(br, dots, 1, 1, 0, nbins, 0, 0, 0, ymax);
