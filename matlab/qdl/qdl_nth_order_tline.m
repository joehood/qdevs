clear variables

dqmin = .1;
dqmax = .1;
dqerr = 0.001;

nnodes = 4; % must be atleast 2
nodes = QdlNode.empty(0);
branches = QdlBranch.empty(0);

sys = QdlSystem(dqmin, dqmax, dqerr);

nodes(1) = QdlNode('node1', 0, 0, 0);

% dc:
if 1
nodes(1).source_type = QdlSystem.SourceDC;
nodes(1).Vdc = 10.0;
end   

% pwm:
if 0
nodes(1).source_type = QdlSystem.SourcePWM;
nodes(1).freq = 10000;
nodes(1).duty = 0.5;
nodes(1).V1 = 0; 
nodes(1).V2 = 1;
end

% sine:
if 0
nodes(1).source_type = QdlSystem.SourceSINE;
nodes(1).freq = 10;
nodes(1).Va = 1;
end

sys.add_node(nodes(1));

C = 1e-6;
L = 1e-6;
s = 2;

for k=2:nnodes-1
    nn = strcat('node', num2str(k));
    nodes(k) = QdlNode(nn, C*10^(k*s), 1, 0);
    nb = strcat('branch', num2str(k-1), '-', num2str(k));
    branches(k-1) = QdlBranch(nb, nodes(k-1), nodes(k), L*10^(k*s), 0.01, 0);
    sys.add_node(nodes(k));
    sys.add_branch(branches(k-1));
end

nodes(nnodes) = QdlNode(strcat('node', num2str(nnodes)), C*10^(nnodes*s), 1, 0);
nb = strcat('branch', num2str(nnodes-1), '-', num2str(nnodes));
branches(nnodes-1) = QdlBranch(nb, nodes(nnodes-1), nodes(nnodes), L*10^(nnodes*s), 1, 0);
sys.add_node(nodes(nnodes));
sys.add_branch(branches(nnodes-1));  

C*10^(nnodes*s)
L*10^(nnodes*s)

sys.init();
%sys.runto(0.1);
sys.runto(1000);
sys.Xdc(1) = 8;
sys.runto(2000);

figure;

for k=1:length(nodes)
    subplot(length(nodes), 2, k*2-1);
    yyaxis left
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k)), 'b-'); hold on;
    %plot(sys.tout(k,1:2:sys.iout(k)), sys.qout(k,1:2:sys.iout(k)), 'k.');
    ylabel('v (V)');
    yyaxis right
    plot(sys.tupd(k,1:sys.iupd(k)), cumsum(sys.nupd(k,1:sys.iupd(k)))/2, 'r-');
    title(nodes(k).name);
    ylabel('updates');
end

xlabel('t (s)');

for k=1:length(branches)
    kb = length(nodes) + k;
    subplot(length(nodes), 2, k*2);
    yyaxis left
    plot(sys.tout(kb,1:sys.iout(kb)), sys.qout(kb,1:sys.iout(kb)), 'b-'); hold on;
    %plot(sys.tout(kb,1:2:sys.iout(kb)), sys.qout(kb,1:2:sys.iout(kb)), 'k.');
    ylabel('i (A)');
    yyaxis right
    plot(sys.tupd(kb,1:sys.iupd(kb)), cumsum(sys.nupd(kb,1:sys.iupd(kb)))/2, 'r-');
    title(branches(k).name);
    ylabel('updates');
    
    
%     kb = length(nodes) + k;
%     subplot(length(nodes), 2, k*2+2);    
%     plot(sys.tout(kb,1:sys.iout(kb)), sys.qout(kb,1:sys.iout(kb)), 'r-'); hold on;
%     plot(sys.tout(kb,1:sys.iout(kb)), sys.qout(kb,1:sys.iout(kb)), 'k.');
%     title(strcat(branches(k).name, "(zoom)"));
%     ylabel('i (A)');
end

xlabel('t (s)');
