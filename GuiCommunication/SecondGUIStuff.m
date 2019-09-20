classdef SecondGUIStuff < handle
    
    events
        StuffChanged
    end
    
    properties
         % Default parameters set so object initialises in a valid state
        parameter1 = 'SomeText';
        parameter2 = 0.57;
    end
    
    methods
        function set.parameter1( obj, param1 )
            validateattributes( param1, { 'char' }, { } )
            
            obj.parameter1 = param1;
        end
        
        function set.parameter2( obj, param2 )
            validateattributes( param2, { 'double' }, { 'scalar', 'nonnegative', '<=', 1 } )
            
            obj.parameter2 = param2;
        end
    end
end