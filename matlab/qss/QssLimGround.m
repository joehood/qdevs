

classdef QssLimGround < QssLimNode 
  
    methods
  
        function self = QssLimGround(v0)
            
            self@QssLimNode('gnd', 0.0, 0.0, 0.0, v0);
            
            self.is_gnd = 1;
            
        end
        
    end
    
end