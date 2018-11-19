

classdef QssOdeAtom < handle

    properties
        x      % internal state
        d      % derivative
        q      % external (quantized) state
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
    end  
  
    methods
  
        function obj = QssOdeAtom(dQ, eps, a, b, c)
            obj.dQ = dQ;
            obj.eps = eps; 
            obj.a = a;
            obj.b = b;
            obj.c = c;
        end
        
        function init(obj)
            obj.x = 0;
            obj.q = 0;
            obj.d = 0;
            obj.qlast = 0;
            obj.tlast = 0;
            obj.tnext = 0;
            obj.thist(1) = 0;
            obj.qhist(1) = 0;
            obj.k = 2;  
        end

        function trigger = update(obj, input, tstop)

            trigger = 0;  % initialize trigger flag to false
            
            % 1. transition (delta_int / delta_ext) function (update x):
            
            t = obj.tnext;      % current time is tnext
            e = t - obj.tlast;  % elapsed time
            obj.tlast = t;      % save tlast for next elapsed time calculation
            
            obj.d = (obj.c - obj.q * obj.b + input) / obj.a;  % static function
            obj.x = obj.x + obj.d * e;                        % integrator
                     
            % 2. output function (lambda) function (update q):
            
            if obj.x >= obj.q + obj.dQ - obj.eps
                obj.q = obj.q + obj.dQ;
            elseif obj.x <= obj.q - 0.5 * obj.dQ + obj.eps
                obj.q = obj.q - obj.dQ;
            end

            if obj.q ~= obj.qlast
            
                obj.qlast = obj.q;  % save qlast for next comparison
                trigger = 1;        % trigger ext event
                
                % save data:
                obj.thist(obj.k) = obj.tlast;
                obj.qhist(obj.k) = obj.q; 
                obj.k = obj.k + 1;
                
            end

            % 4. time advance (ta) function:
            
            if (obj.d > 0)
                obj.tnext = t + (obj.q + obj.dQ - obj.x) / obj.d;
            elseif (obj.d < 0)
                obj.tnext = t + (obj.q - 0.5 * obj.dQ - obj.x) / obj.d;
            else
                obj.tnext = inf;
            end
            
            %obj.tnext = max(obj.tnext, 0);
            obj.tnext = min(obj.tnext, tstop);
            
        end
        
    end
    
end