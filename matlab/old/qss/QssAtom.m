

classdef QssAtom < handle

    properties
        
        name;    % node name (helpful for debugging)
        x;       % internal state
        x0;      % internal state at system time = time0
        d;       % derivative
        q;       % external (quantized) state
        qlast;   % ext state at last internal transition
        dQ;      % quantum
        epsilon; % hysteresis width
        t0;      % initial time
        tlast;   % time at last internal transition
        tnext;   % next transition time (from ta or ext trigger)
        thist;   % time history array
        qhist;   % q history array
        khist;   % history array index
        time;    % simulation time
        trigger; % external update trigger
        
    end  
  
    methods
  
        function self = QssAtom(name, x0)
            
            self.name = name;
            self.x0 = x0;
            
        end
        
        function init(self, t0)
            
            self.t0 = t0;
            self.x = self.x0;
            self.q = self.x0;
            self.qlast = self.x0;
            self.time = t0;
            self.tlast = t0;
            self.d = 0.0;
            self.tnext = inf;            
            self.khist = 1;
            self.thist(1) = self.time; 
            self.qhist(1) = self.q;       
            
        end
        
        function update(self, time)   
            
            self.time = time;            
            self.trigger = 0;

            self.dint();
            self.quantize();
            self.update_derivative(); 
            self.ta();
            
            if self.q ~= self.qlast 
                self.save_history();
                self.dext();
            end
 
        end
        
        function quantize(self)  
            
            self.qlast = self.q;
            
            if self.x >= self.q + self.dQ - self.epsilon
                self.q = self.q + self.dQ;
            elseif self.x <= self.q - 0.5 * self.dQ + self.epsilon
                self.q = self.q - self.dQ;
            end
            
        end
        
        function dint(self)
            
            self.x = self.x + self.d * (self.time - self.tlast);
            self.tlast = self.time;
            
        end

        function ta(self)

           dt = 1.0e-15;
           
           if (self.d > 0.0)
               dt = (self.q + self.dQ - self.x) / self.d;
           elseif (self.d < 0.0)
               dt = (self.q - 0.5 * self.dQ - self.x) / self.d;
           else
               dt = inf;
           end
           
           self.tnext = self.time + abs(dt);
            
        end
        
        function save_history(self) 
            
            if self.thist(self.khist) ~= self.time  % avoid dup. time points
                self.khist = self.khist + 1;
                self.thist(self.khist) = self.time;
            end
                
            self.qhist(self.khist) = self.q; 
        end
        
    end  % end methods
    
    methods (Abstract)

        update_derivative(self);
        
        dext(self);

    end
    
end