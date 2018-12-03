
dqmin = 0.01;
dqmax = 0.01;
dqerr = 0.01;

R = [ 1,    1,    1,    1   ];
L = [ 1e-6, 1e-2, 1e2,  1e6 ];
C = [ 1e-6, 1e-2, 1e2,  1e6 ];
G = [ 1,    1,    1,    1   ];

sys = QdlSystem(dqmin, dqmax, dqerr);

n = 4; % must be even

nodes = QdlNode.empty(0);
branches = QdlBranch.empty(0);

nn = 0;
nb = 0;

for i=1:n
    
    for j=1:n
        
        if i <= n/2
            if j <= n/2
                r = R(1); l = L(1); c = C(1); g = G(1);
            else
                r = R(3); l = L(3); c = C(3); g = G(3);
            end
        else
            if j <= n/2
                r = R(2); l = L(2); c = C(2); g = G(3);
            else
                r = R(4); l = L(4); c = C(4); g = G(4);
            end
        end

        nn = (i-1)*n+j;
        lbl = strcat('node', num2str(nn), '(C=', num2str(l), ')');
        disp(num2str(nn));
        nodes(nn) = QdlNode(lbl, c, g, 0);
        sys.add_node(nodes(nn));
        
        l = 1;
        r = 0.01;
        
        if i > 1
            nb = nb + 1;
            lbl = strcat('branch', num2str(nb), '(L=', num2str(l), ')');
            branches(nb) = QdlBranch(lbl, l, r, 0);
            ni = (i-2)*n+j;
            nj = nn;
            disp(strcat(num2str(ni),'-',num2str(nj)));
            branches(nb).connect(nodes(ni), nodes(nj));
            sys.add_branch(branches(nb));
            
        end
        
        if j > 1
            nb = nb + 1;
            lbl = strcat('branch', num2str(nb), '(L=', num2str(l), ')');
            branches(nb) = QdlBranch(lbl, l, r, 0);
            ni = (i-1)*n+(j-1);
            nj = nn;
            disp(strcat(num2str(ni),'-',num2str(nj)));
            branches(nb).connect(nodes(ni), nodes(nj));
            sys.add_branch(branches(nb));
            
        end
        
    end
    
end

sys.init();
sys.runto(100);
sys.H(1) = 10.0;
sys.runto(15000);

nbins = 500;
ymax = 0;

disp(strcat('order=', num2str(nn+nb)));

figure;

subplot(4, 1, 1);
sys.plot(nodes(1), 0, 1, 1, 0, nbins, 0, 0, 0, ymax);

subplot(4, 1, 2);
sys.plot(nodes(1+n), 0, 1, 1, 0, nbins, 0, 0, 0, ymax);

subplot(4, 1, 3);
sys.plot(nodes(nn-n), 0, 1, 1, 0, nbins, 0, 0, 0, ymax);

subplot(4, 1, 4);
sys.plot(nodes(nn), 0, 1, 1, 0, nbins, 0, 0, 0, ymax);



