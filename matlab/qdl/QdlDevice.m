

classdef QdlDevice < handle

    properties

        name
        dqmax
        dqmin
        dqerr     
        index
        freq
        duty
        phi

    end

    methods

        function self = QdlDevice(name)

            self.name  = name;
            self.dqmax = 0;
            self.dqmin = 0;
            self.dqerr = 0;     
            self.index = 0;
            self.freq  = 0;
            self.duty  = 0;
            self.phi   = 0;

        end

    end

end
