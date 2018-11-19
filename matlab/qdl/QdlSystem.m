


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