classdef kmpc_yalmip < matlab.System & matlab.system.mixin.Propagates
    % An MPC controller implemented using YALMIP
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        Ts = 2e-3
        ny = 3
        N = 5
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)
        controller
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
%             loadStruct = load('data/kmpc_yalmip.mat','controller');
%             obj.controller = loadStruct.controller;
            load('data/kmpc_yalmip.mat','mpc_setup');
            obj.controller = kmpc_setup_y2f(mpc_setup);
        end

        function [uopt,yopt,exitflag,solvetime] = stepImpl(obj,z0,T_ref)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            [solution,exitflag,~,~,~,info] = obj.controller({z0,T_ref});
            uopt = solution{1};
            yopt = solution{2};
            solvetime = info.solvertime;
        end
        
        function sts = getSampleTimeImpl(obj)
            % Set the sample time
            sts = createSampleTime(obj,'Type','Discrete',...
                                       'SampleTime',obj.Ts,...
                                       'OffsetTime',0);
            %sts = createSampleTime(obj,'Type','Fixed In Minor Step');
        end
        
        function resetImpl(~)
            % Initialize / reset discrete-state properties
        end
        
        function num = getNumInputsImpl(~)
            num = 2;
        end
        function num = getNumOutputsImpl(~)
            num = 4;
        end
        function [dt1,dt2,dt3,dt4] = getOutputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
            dt4 = 'double';
        end
        function [dt1,dt2] = getInputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
        end
        function [sz1,sz2,sz3,sz4] = getOutputSizeImpl(obj)
            sz1 = [1,obj.N];        % uopt
            sz2 = [obj.ny*obj.N,1]; % predictions
            sz3 = [1,1];            % exitflag
            sz4 = [1,1];            % solve time
        end
        function sz1 = getInputSizeImpl(~)
        	sz1 = [1,1];
        end
        function cp1 = isInputComplexImpl(~)
        	cp1 = false;
        end
        function [cp1,cp2,cp3,cp4] = isOutputComplexImpl(~)
        	cp1 = false;
            cp2 = false;
            cp3 = false;
            cp4 = false;
        end
        function fz1 = isInputFixedSizeImpl(~)
        	fz1 = true;
        end
        function [fz1,fz2,fz3,fz4] = isOutputFixedSizeImpl(~)
        	fz1 = true;
            fz2 = true;
            fz3 = true;
            fz4 = true;
        end
    end
end
