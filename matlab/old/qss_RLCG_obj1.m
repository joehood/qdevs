clear all

R = 1;
L = 1;
E = 1;
G = 1;
C = 1;
H = 1;

t = 0.0;
tstop = 10.0;

dQ = 0.01;
eps = dQ/4;
branch = QssAtom(dQ, eps);
node = QssAtom(dQ, eps);
branch.init();
node.init();
branch.d = fbranch(R, L, E, branch.q, node.q);
node.d = fnode(G, C, H, branch.q, node.q);
t = branch.update(0);
t = node.update(0);

while (t < tstop)

    %disp(t);
    
    if branch.tnext == node.tnext
    
        i = branch.q;
        v = node.q;
    
        branch.d = fbranch(R, L, E, i, v);
        
        [t, trigger] = branch.update(t); 
        if trigger
            node.tnext = t;
        end

        node.d = fnode(G, C, H, i, v);
        
        [t, trigger] = node.update(t);
        if trigger 
            branch.tnext = t;
        end
        
    elseif branch.tnext < node.tnext 

        branch.d = fbranch(R, L, E, branch.q, node.q); 
        
        [t, trigger] = branch.update(t);
        if trigger
            node.tnext = t;
        end
        
    elseif branch.tnext > node.tnext 
    
        node.d = fnode(G, C, H, branch.q, node.q);
        if node.update() 
            branch.tnext = t;
        end
        
    end
    
end


disp(branch.k)

plot(branch.thist(1:branch.k-2), branch.qhist(1:branch.k-2), 'b-o'); hold on;
plot(node.thist(1:node.k-1), node.qhist(1:node.k-1), 'r-o');

%disp('done.')

function dbranch = fbranch(R, L, E, qbranch, qnode)
    dbranch = 1/L * (E - qbranch*R - qnode);
end

function dnode = fnode(G, C, H, qbranch, qnode)
    dnode = 1/C * (H - qnode*G + qbranch);
end


