

classdef QssAtom < handle

    properties
        x;      % internal state
        x0;     % internal state at system time = time0
        d;      % derivative
        q;      % external (quantized) state
        qlast;  % ext state at last internal transition
        dQ;     % quantum
        eps;    % hysteresis width
        tlast;  % time at last internal transition
        tnext;  % next transition time (from ta or ext trigger)
        thist;  % time history array
        qhist;  % q history array
        k;      % history array index
    end  
  
    methods
  
        function obj = QssAtom(dQ, eps)
            obj.dQ = dQ;
            obj.eps = eps; 
        end
        
        function init(obj, t0)
            
            obj.x = obj.x0;
            obj.q = obj.x0;
            obj.qlast = obj.x0;
            obj.tlast = t0;
            
            obj.d = 0;
            obj.tnext = inf;
            
            obj.thist(1) = t0;
            obj.qhist(1) = obj.q;
            obj.k = 1; 
            
            obj.update_derivative();
            
        end
        
        function quantize(obj)  
            
            if obj.x >= obj.q + obj.dQ - obj.eps
                obj.q = obj.q + obj.dQ;
            elseif obj.x <= obj.q - 0.5 * obj.dQ + obj.eps
                obj.q = obj.q - obj.dQ;
            end
            
        end

        function advance(obj, t)

           if (obj.d > 0)
               dt = (obj.q + obj.dQ - obj.x) / obj.d;
               obj.tnext = t + abs(dt);
           elseif (obj.d < 0)
               dt = (obj.q - 0.5 * obj.dQ - obj.x) / obj.d;
               obj.tnext = t + abs(dt);
           else
               obj.tnext = inf;
           end
            
        end
        
        
        function transition(obj, t)

            dt = t - obj.tlast;
            obj.tnext = t;
            obj.x = obj.x + obj.d * dt;

        end
        
    end  % end methods
    
    methods (Abstract)

        update_derivative(obj);

    end
    
end