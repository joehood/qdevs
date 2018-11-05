

classdef QssLimNode < QssAtom

    properties
        
        H;        % shunt ideal current
        G;        % shunt conductance
        C;        % shunt capacitance
        invC;     % 1/C
        is_gnd;   % is ground/reference node
        branches; % handles to connected branches
        nbranch;  % number of connected branches
        brn_sign; % branch direction signs
        
    end  
  
    methods
  
        function obj = QssLimNode(H, G, C, v0, is_gnd, dQ, eps)
            
            obj@QssAtom(dQ, eps);
            
            obj.H = H;
            obj.G = G;
            obj.C = C;
            obj.invC = 1/C;
            obj.is_gnd = is_gnd;
            
            obj.x0 = v0;     
            
            obj.nbranch = 0;
            
            obj.branches = QssLimBranch.empty(0);
            obj.brn_sign = double.empty(0);
        end
        
        function connect_branch(obj, branch, direction)
            
            obj.nbranch = obj.nbranch + 1;
            obj.branches(obj.nbranch) = branch;
            obj.brn_sign(obj.nbranch) = direction;
         
        end
        
        function update_derivative(obj)
            
            isum = sum(obj.branches.q * obj.brn_sign);
            
            obj.d = obj.invC * (obj.H - obj.x * obj.G - isum);
        
        end
        
    end
    
end