

classdef LiqssLimGround < LiqssLimNode 
  
    methods
  
        function self = LiqssLimGround(name, v0)
            
            self@LiqssLimNode(name, 0.0, 0.0, 0.0, v0, 0.0, 0.0);
            
            self.constant_voltage = 1;
            
        end
        
    end
    
end