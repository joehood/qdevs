clear variables

dq = 0.001;

sys = QdlSystem(dq, dq, 0.001);

vg = QdlNode('vg', 0, 0, 0);

vg.source_type = QdlSystem.SourcePWM;
vg.v1 = 0.0;
vg.v2 = 1.0;
vg.freq = 10;
vg.duty = 0.5;

v1 = QdlNode('v1', 0.1, .1, 0);

b1 = QdlBranch('b1', 1, 1, 0);

b1.connect(vg, v1);

sys.add_node(vg);
sys.add_node(v1);
sys.add_branch(b1);

sys.init();
sys.runto(1);
sys.xdc(1) = 1;
sys.runto(10);

figure;

subplot(3, 1, 1);

for k=2:sys.n
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k))); hold on;
end

ylabel('atom x');
title(strcat('dq=', num2str(dq)));

dq = 0.001;

v1.dqmin = dq;
v1.dqmax = dq;
b1.dqmin = dq;
b1.dqmax = dq;

sys.init();
sys.runto(1);
sys.xdc(1) = 1;
sys.runto(10);

subplot(3, 1, 2);

for k=2:sys.n
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k))); hold on;
end
ylabel('atom x');
title(strcat('dq=', num2str(dq)));

dq = 0.05;

v1.dqmin = dq;
v1.dqmax = dq;
b1.dqmin = dq;
b1.dqmax = dq;

sys.init();
sys.runto(1);
sys.xdc(1) = 1;
sys.runto(10);

subplot(3, 1, 3);

for k=2:sys.n
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k))); hold on;
end
ylabel('atom x');
title(strcat('dq=', num2str(dq)));

xlabel('t (s)');


