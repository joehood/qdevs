

classdef LiqssLimNode < LiqssAtom

    properties
      
        H;        % shunt ideal current
        G;        % shunt conductance
        C;        % shunt capacitance
        den;      % 1/C
        branches; % handles to connected branches
        nbranch;  % number of connected branches
        brn_sign; % branch direction signs
        
        constant_voltage;   % is ground/reference node (x = q = x0)
        
    end  
  
    methods
  
        function self = LiqssLimNode(name, H, G, C, v0, dQVoltage, epsVoltage)
            
            self@LiqssAtom(name, v0, dQVoltage, epsVoltage);
            
            self.H = H;
            self.G = G;
            self.C = C;
            self.den = 1/C;
            self.constant_voltage = 0;              
            self.nbranch = 0;         
            self.branches = LiqssLimBranch.empty(0);
            self.brn_sign = double.empty(0);
            
        end
        
        function connect_branch(self, branch, direction)
            
            self.nbranch = self.nbranch + 1;
            self.branches(self.nbranch) = branch;
            self.brn_sign(self.nbranch) = direction;
         
        end
        
        function [d] = f(self, q)
            
            isum = 0;
            
            for k = 1:self.nbranch
                isum = isum + self.branches(k).q * self.brn_sign(k);
            end
            
            d = self.den * (self.H - q * self.G - isum);
        
        end
       
        function dext(self)
           
           for k = 1:self.nbranch                 
               self.branches(k).trigger = 1;           
           end 
           
        end
        
    end
    
end