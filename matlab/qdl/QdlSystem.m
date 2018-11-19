
%
%  n   := number of nodes
%  m   := number of branches
%
%  a(n, m) := incidence matrix between node i and branch k
%             0: no connection
%             1: ith node of branch k connected to node i
%            -1: jth node of branch k connected to node i
%
%  b(n, n) := vccs gain matrix between node i and node p
%  s(n, m) := cccs gain matrix between node i and branch k
%  t(m, n) := vcvs gain matrix between branch k and node i
%  z(m, m) := ccvs gain matrix between branch k and branch q 
%
%  c(n) := capacitance at node i
%  g(n) := conductance at node i
%  h(n) := current injection at node i
%
%  l(n) := series inductance in branch k
%  r(n) := series resistance in branch k
%  e(n) := voltage at branch k from node i to node j
%
%  where:
% 
%  vi = 1/ci * (hi + bip * vp + siq * iq + sum(ik) - gi * vi)
%
%  ik = 1/lk * (ek + tkp * vp + zkq * iq + (vik - vjk) - rk * ik)
%
%


classdef QdlSystem < handle
    
    properties
        
        an;   % node atoms
        ab;   % branch atoms
        nn;   % number of atoms (size of state arrays)
        m;    % num history values for state arrays
        tmax; % max simulation time
        t;    % current simulation time
        x;    % internal states
        q;    % quantum states
        d;    % internal state derivative
        dQ;   % quantum levels

    end
    
    methods
        
        function self = QdlSystem()
            
            self.n = 0;
            self.a = QdlAtom.empty(0);
            self.m = 2;
            
        end
        
        function add(self, a)
            
            self.n = self.n + 1;
            self.a(self.n) = a;
            
        end
        
        function init(self, tmax)
            
            self.tmax = tmax;
            self.x = zeros(self.n, self.m);
            self.q = zeros(self.n, self.m);
            self.d = zeros(self.n, self.m);
            
        end
        
    end

end