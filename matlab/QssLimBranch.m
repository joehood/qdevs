

classdef QssLimBranch < QssAtom

    properties
        
        nodei; % node handle for node at terminal i
        nodej; % node handle for node at terminal j
        E;     % series ideal voltage
        R;     % series resistance
        L;     % series inductance
        invL;  % 1/L
        
    end  
  
    methods
  
        function obj = QssLimBranch(E, R, L, i0, nodei, nodej, dQ, eps)
            
            obj@QssAtom(dQ, eps);
            
            obj.nodei = nodei;
            obj.nodej = nodej;
            obj.E = E;
            obj.R = R;
            obj.L = L;
            obj.invL = 1/L;
            
            obj.x0 = i0;
            
            obj.nodei = nodei;
            obj.nodej = nodej;
            
            obj.nodei.connect_branch(obj, 1);
            obj.nodej.connect_branch(obj, -1);
            
        end
        
        function update_derivative(obj)
            
            obj.d = obj.invL * (obj.E - obj.x * obj.R - obj.nodei.x - obj.nodej.x);
        
        end
        
    end
    
end