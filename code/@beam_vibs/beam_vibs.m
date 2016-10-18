classdef beam_vibs < handle
%% classdef beam_vibs
% 
% 
% 
% author: john devitis
% create date: 18-Oct-2016 16:46:00

%% object properties
	properties
        K % global stiffness matrix
        M % global mass matrix
        V % eigenvectors
        W % natural frequencies [rad/sec]
        F % natural frequencies [hz]
        Mr % modal mass
        Kr % modal stiffness
        Vn % mass normalized shapes
        dampr % damping ratio
        ns % number of spectral lines
        freqbnds % frequency bnds
        outLoc 
        inLoc
        AA
        HH
	end

%% dependent properties
	properties (Dependent)
        no
        ni
        ne % number of effective modes
        w
        dampf % damping factor [rad/sec]
        Wn % damped natural frequency [rad/sec]
        root % positive system poles
        Qr
	end

%% private properties
	properties (Access = private)
	end

%% dynamic methods
	methods
	%% constructor
		function self = beam_vibs()
		end

	%% ordinary methods

	%% dependent methods
        function ne = get.ne(self)
            ne = size(self.V,2);
        end
	end

%% static methods
	methods (Static)
	end

%% protected methods
	methods (Access = protected)
	end

end
