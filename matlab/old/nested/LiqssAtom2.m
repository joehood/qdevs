

classdef LiqssAtom2 < handle

    properties
        
        name;    % node name (helpful for debugging)
        x;       % internal state
        x0;      % internal state at system time = time0
        d;       % derivative
        dlast;   % last derivative
        q;       % external (quantized) state
        qhi;     % upper quantized state
        qlo;     % lower quantized state
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
        dtmin;   % minimum time advance
        
    end  
  
    methods
  
        function self = LiqssAtom2(name, x0, dQ, epsilon)
            
            self.name = name;
            self.x0 = x0;
            self.dQ = dQ
            self.epsilon = epsilon;
            self.dtmin = 1e-12;
            
        end
        
        function init(self, t0)
            
            self.t0 = t0;
            self.x = self.x0;
            self.q = self.x0;
            self.qlast = self.x0;
            self.qhi = self.q + self.dQ;
            self.qlo = self.q - self.dQ;
            self.time = t0;
            self.tlast = t0;
            self.d = 0.0;
            self.dlast = 0.0;
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
            self.d = self.f(self.q);
            self.ta();
            
            if self.q ~= self.qlast 
                self.qlast = self.q;
                self.save_history();
                self.dext();
            end
 
        end
        
        function quantize(self)  
            
            %qhi = self.qhi;
            %qlo = self.qlo;
            %qlast = self.q;
            
            self.dlast = self.d;
            
            change = 0;
            
            if self.x >= self.qhi
                self.q = self.qhi;
                self.qlo = self.qlo + self.dQ;
                change = 1;
            elseif self.x <= self.qlo
                self.q = self.qlo;
                self.qlo = self.qlo - self.dQ;
                change = 1;
            end
            
            self.qhi = self.qlo + 2 * self.dQ;
            
            if change  % we've ventured out of qlo/qhi bounds
            
                self.d = self.f(self.q);
                
                % if the derivative has changed signs, then we know 
                % we are in a potential oscillating situation, so
                % we will set the q such that the derivative=0:
                
                if self.d * self.dlast < 0  % if derivative has changed sign
                    flo = self.f(self.qlo); 
                    fhi = self.f(self.qhi);
                    k = (2 * self.dQ) / (fhi - flo);
                    self.q = self.qhi - k * fhi;
                end
                
            end
            
        end
        
        function dint(self)
            
            self.x = self.x + self.d * (self.time - self.tlast);
            self.tlast = self.time;
            
        end

        function ta(self)

            % estimate the time to the next qhi/qlo crossing:
            if (self.d > 0)
                self.tnext = self.time + (self.qhi - self.x) / self.d;
            elseif (self.d < 0)
                self.tnext = self.time + (self.qlo - self.x) / self.d;
            else
                self.tnext = inf;
            end
            
            % tnext from internal event shall not be 0:
            self.tnext = max(self.tnext, self.tlast + self.dtmin);
            
            % also force an update at tstop if tnext > tstop:
            %self.tnext = min(self.tnext, self.tstop);
            
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

        d = f(self, q);
        
        dext(self);

    end
    
end