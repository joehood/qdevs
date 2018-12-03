

classdef QssLimBranch < QssAtom

    properties
        
        nodei; % node handle for node at terminal i
        nodej; % node handle for node at terminal j
        E;     % series ideal voltage
        R;     % series resistance
        L;     % series inductance
        den;   % 1/L
        
    end  
  
    methods
  
        function self = QssLimBranch(name, E, R, L, i0, nodei, nodej)
            
            self@QssAtom(name, i0);
            
            self.nodei = nodei;
            self.nodej = nodej;
            self.E = E;
            self.R = R;
            self.L = L;
            self.den = 1/L;           
            self.nodei = nodei;
            self.nodej = nodej;          
            self.nodei.connect_branch(self, 1);
            self.nodej.connect_branch(self, -1);
            
        end
        
        function update_derivative(self)
            
            if self.nodei.is_gnd && self.nodej.is_gnd
                vij = 0.0;
            elseif self.nodei.is_gnd
                vij = -self.nodej.q;
            elseif self.nodej.is_gnd
                vij = self.nodei.q;
            else
                vij = self.nodei.q - self.nodej.q;
            end
            
            self.d = self.den * (self.E - self.q * self.R + vij);
        
        end
       
        function dext(self)
           
           self.nodei.trigger = 1;  
           self.nodej.trigger = 1;

        end 
        
    end
    
end