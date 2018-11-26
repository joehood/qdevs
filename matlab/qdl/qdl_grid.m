clear variables

dqmin = 0.00001;
dqmax = 0.001;
maxerr = 0.01;

npt = 1e5;

%sys.add_node(i, C, G, H, p, B, q, S, dqmin, dqmax, maxerr);
%sys.add_branch(k, i, j, L, R, E, p, T, q, Z, dqmin, dqmax, maxerr);

n = 10;

% source branch:

sys = QdlSystem(n*n, n*n*2, npt);

L = 0.1;
R = 1.0;
C = 0.1;
G = 1.0;
H = 1.0;

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
        
        sys.add_node(i*j, C, G, 0.0, 0, 0.0, 0, 0.0, dqmin, dqmax, maxerr);      
        sys.add_branch(i*j, i+j, i+j-1, L, R, 0.0, 0, 0.0, 0, 0.0, dqmin, dqmax, maxerr);
        sys.add_branch(i*j*2, i+j-1, i+j, L, R, 0.0, 0, 0.0, 0, 0.0, dqmin, dqmax, maxerr);

    end
 
end

sys.init();
sys.runto(100);
sys.E(1) = 1.0;
sys.runto(5000);
sys.E(1) = 2.0;
sys.runto(10000);

for k = 1:10:n*n+n*n*2
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k)), 'b.-'); hold on;
end

figure;
for k = 1:10:n*n+n*n*2
    plot(sys.tupd(k, 1:sys.iupd(k)), cumsum(sys.nupd(k, 1:sys.iupd(k)))); hold on;
end

