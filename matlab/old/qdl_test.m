clear all

R = 1;
L = 1;
E = 1;
G = 1;
C = 1;
H = 1;

t = 0.0;
tstop = 10.0;

dQ = 0.001;
eps = dQ/4;

branch = QssOdeAtom(dQ, eps, L, R, E);
node = QssOdeAtom(dQ, eps, C, G, H);

branch.init();
node.init();

t = branch.update(0, tstop);
t = node.update(0, tstop);

while (t < tstop)

    isum = branch.q;
    vij = -node.q;  % (0 - v_node)
        
    if branch.tnext == node.tnext
    
        [t, trigger] = branch.update(vij, tstop); 
        if trigger node.tnext = t; end

        [t, trigger] = node.update(isum, tstop);
        if trigger branch.tnext = t; end
        
    elseif branch.tnext < node.tnext 

        [t, trigger] = branch.update(vij, tstop);
        if trigger node.tnext = t;  end
        
    elseif branch.tnext > node.tnext 

        [t, trigger] = node.update(isum, tstop);
        if trigger branch.tnext = t; end
        
    end
    
end


disp(branch.k)

plot(branch.thist, branch.qhist, 'b-o'); hold on;
plot(node.thist, node.qhist, 'r-o');

%disp('done.')

