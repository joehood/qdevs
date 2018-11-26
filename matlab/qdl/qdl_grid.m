clear variables

dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.01;

sys = QdlSystem(dqmin, dqmax, dqerr);

L = 0.1;
R = 1.0;
C = 0.1;
G = 1.0;
H = 1.0;

nodes = QdlNode.empty(0);
branches = QdlBranch.empty(0);
kn = 0;
kb = 0;

n = 4;

kn = kn + 1;
nodes(kn) = QdlNode('', C, G, 0.0);
sys.add_node(nodes(kn));

for i = 1:n

    for j = 1:n
        
        if i < n/2
            if j < n/2
                C = 1000.0;
                L = 1000.0;
            else
                C = 0.1;
                L = 0.1;
            end  
        else
            if j < n/2
                C = 10.0;
                L = 10.0;
            else
                C = 0.001;
                L = 0.001;
            end
        end
        
        kn = kn + 1;
        nodes(kn) = QdlNode('', C, G, 0.0);
        sys.add_node(nodes(kn));
        
        kb = kb + 1;
        branches(kb) = QdlBranch('', nodes(i+j), nodes(i+j-1), L, R, 0.0);
        sys.add_branch(branches(kb));
        
        kb = kb + 1;
        branches(kb) = QdlBranch('', nodes(i+j-1), nodes(i+j), L, R, 0.0);
        sys.add_branch(branches(kb));

    end
 
end

sys.init();
sys.runto(100);
sys.E(1) = 1.0;
sys.runto(500);
sys.E(1) = 2.0;
sys.runto(1000);

for k = 1:10:n*n+n*n*2
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k)), 'b.-'); hold on;
end

figure;

for k = 1:10:n*n+n*n*2
    plot(sys.tupd(k, 1:sys.iupd(k)), cumsum(sys.nupd(k, 1:sys.iupd(k)))); hold on;
end

