

classdef LiqssLimBranch < LiqssAtom

    properties
        
        nodei; % node handle for node at terminal i
        nodej; % node handle for node at terminal j
        E;     % series ideal voltage
        R;     % series resistance
        L;     % series inductance
        den;   % 1/L
        
    end  
  
    methods
  
        function self = LiqssLimBranch(name, E, R, L, i0, nodei, nodej, dQCurrent, epsCurrent)
            
            self@LiqssAtom(name, i0, dQCurrent, epsCurrent);
            
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
        
        function [d] = f(self, q)
            
            if self.nodei.constant_voltage && self.nodej.constant_voltage
                vij = 0.0;
            elseif self.nodei.constant_voltage
                vij = -self.nodej.q;
            elseif self.nodej.constant_voltage
                vij = self.nodei.q;
            else
                vij = self.nodei.q - self.nodej.q;
            end
            
            d = self.den * (self.E - q * self.R + vij);
        
        end
       
        function dext(self)
           
           self.nodei.trigger = 1;  
           self.nodej.trigger = 1;

        end 
        
    end
    
end