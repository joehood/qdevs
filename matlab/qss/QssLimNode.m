

classdef QssLimNode < QssAtom

    properties
      
        H;        % shunt ideal current
        G;        % shunt conductance
        C;        % shunt capacitance
        den;      % 1/C
        is_gnd;   % is ground/reference node
        branches; % handles to connected branches
        nbranch;  % number of connected branches
        brn_sign; % branch direction signs
        
    end  
  
    methods
  
        function self = QssLimNode(name, H, G, C, v0)
            
            self@QssAtom(name, v0);
            
            self.H = H;
            self.G = G;
            self.C = C;
            self.den = 1/C;
            self.is_gnd = 0;              
            self.nbranch = 0;         
            self.branches = QssLimBranch.empty(0);
            self.brn_sign = double.empty(0);
            
        end
        
        function connect_branch(self, branch, direction)
            
            self.nbranch = self.nbranch + 1;
            self.branches(self.nbranch) = branch;
            self.brn_sign(self.nbranch) = direction;
         
        end
        
        function update_derivative(self)
            
            isum = 0;
            
            for k = 1:self.nbranch
                isum = isum + self.branches(k).q * self.brn_sign(k);
            end
            
            self.d = self.den * (self.H - self.q * self.G - isum);
        
        end
       
        function dext(self)
           
           for k = 1:self.nbranch                 
               self.branches(k).trigger = 1;           
           end 
           
        end
        
    end
    
end