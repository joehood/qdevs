

classdef LiqssLimGround2 < LiqssLimDevice 
  
    methods
  
        function self = LiqssLimGround2(name, v0)
            
            self@LiqssLimDevice(name, 0.0, 0.0, 0.0);
            
            ground_node = LiqssLimNode2('gnd', 0.0, 0.0, 0.0, v0, 0.0, 0.0)
            
            ground_node.is_gnd = 1;
            
            self.add_node(ground_node);
            
        end
        
    end
    
end