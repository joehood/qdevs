


classdef LiqssOdeAtom < handle

    properties
        x      % internal state
        d      % derivative
        dlast  % last derivative
        q      % external (quantized) state
        qhi    % upper quantized state
        qlo    % lower quantized state
        qlast  % ext state at last internal transition
        dQ     % quantum
        eps    % hysteresis width
        tlast  % time at last internal transition
        tnext  % next transition time (from ta or ext trigger)
        thist  % time history array
        qhist  % q history array
        k      % history array index
        a      % ODE constant a  (where a*x' = b*x + c)
        b      % ODE constant b
        c      % ODE constant c
        MIN_STEP
    end  
  
    methods
  
        function self = LiqssOdeAtom(dQ, eps, a, b, c)
            self.dQ = dQ;
            self.eps = eps; 
            self.a = a;
            self.b = b;
            self.c = c;
            self.MIN_STEP = 1e-15;
        end
        
        function init(self)
            self.x = 0;
            self.d = 0;
            self.dlast = 0; 
            self.q = 0;
            self.qhi = self.q + self.dQ;
            self.qlo = self.q - self.dQ;
            self.qlast = 0;
            self.tlast = 0;
            self.tnext = 0;
            self.thist(1) = 0;
            self.qhist(1) = 0;
            self.k = 2;  
        end

        function trigger = update(self, input, tstop)

            trigger = 0;  % initialize trigger flag to false
            
            % 1. transition (delta_int / delta_ext) function (update x):
            
            t = self.tnext;      % current time is tnext
            e = t - self.tlast;  % elapsed time
            self.tlast = t;      % save tlast for next elapsed time calculation
            
            self.d = (self.c - self.q * self.b + input) / self.a; % (no effect for int transition)
            self.x = self.x + self.d * e; % integrator
                     
            % 2. output function (lambda) function (update q):
            
            qhi = self.qhi;
            qlo = self.qlo;
            qlast = self.q;
            
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
            
                self.d = (self.c - self.q * self.b + input) / self.a;
                
                % if the derivative has changed signs, then we know 
                % we are in a potential oscillating situation, so
                % we will set the q such that the derivative=0:
                
                if self.d * self.dlast < 0  % if derivative has changed sign
                    flo = (self.c - qlo * self.b + input) / self.a; 
                    fhi = (self.c - qhi * self.b + input) / self.a;
                    k = (2 * self.dQ) / (fhi - flo);
                    self.q = qhi - k * fhi;
                end
            end
            
            if self.q ~= self.qlast 
            
                self.qlast = self.q;  % save qlast for next comparison
                trigger = 1;          % trigger ext event
            
                % save data:
                self.thist(self.k) = self.tlast;
                self.qhist(self.k) = self.q; 
                self.k = self.k + 1;
            end

            % 3. time advance (ta) function:
            
            % get the latest estimate of the derivative:
            self.d = (self.c - self.q * self.b + input) / self.a;
            
            % estimate the time to the next qhi/qlo crossing:
            if (self.d > 0)
                self.tnext = t + (self.qhi - self.x) / self.d;
            elseif (self.d < 0)
                self.tnext = t + (self.qlo - self.x) / self.d;
            else
                self.tnext = inf;
            end
            
            % tnext from internal event shall not be 0:
            self.tnext = max(self.tnext, self.tlast+self.MIN_STEP);
            
            % also force an update at tstop if tnext > tstop:
            self.tnext = min(self.tnext, tstop);
            
        end
        
    end
    
end