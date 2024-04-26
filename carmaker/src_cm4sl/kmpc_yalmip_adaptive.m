classdef kmpc_yalmip_adaptive < matlab.System & matlab.system.mixin.Propagates
    % An adaptive MPC controller implemented using YALMIP

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
        controllers
        kmpc_datas
        speed_limits
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            s = load('../../data/kmpc_yalmip.mat','mpc_setup');
            obj.controllers = kmpc_setup_y2f_adaptive(s.mpc_setup);
            s = load('../../models/kmpc_data_adaptive','kmpc_datas','vx');
            obj.kmpc_datas = s.kmpc_datas;
            obj.speed_limits = s.vx;
        end

        function [uopt,yopt,exitflag,solvetime] = stepImpl(obj,vx,ocp_state,kappa_ref,T_ref)
            % find which controller to use
            ctrl_id = find(obj.speed_limits > vx*3.6, 1);
            
            % use the last one if the speed is too high
            if isempty(ctrl_id); ctrl_id = length(obj.speed_limits); end
            
            % normalize the ocp state (not e_int)
            ocp_state_norm = mapstd_custom('apply',ocp_state(1:2),obj.kmpc_datas{ctrl_id}.PX);
            
            % lift the normalized state
            lifted_state = lifting_function(ocp_state_norm,obj.kmpc_datas{ctrl_id});
            
            % add e_int,kappa_ref
            z0 = [lifted_state; ocp_state(3); kappa_ref];
            
            % normalize T_ref
            T_ref = mapstd_custom('apply',T_ref,obj.kmpc_datas{ctrl_id}.PU);
            
            % get the optimal input
            controller = obj.controllers{ctrl_id};
            [solution,exitflag,~,~,~,info] = controller({z0,T_ref});
            
            % set the block outputs
            uopt = solution{1};
            yopt = solution{2};
            solvetime = info.solvertime;
            
            % unnormalize the torque outputs
            uopt = mapstd_custom('reverse',uopt,obj.kmpc_datas{ctrl_id}.PU);
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
            num = 4;
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
        function [dt1,dt2,dt3,dt4] = getInputDataTypeImpl(~)
        	dt1 = 'double';
            dt2 = 'double';
            dt3 = 'double';
            dt4 = 'double';
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
