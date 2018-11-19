classdef PicoFlexCamera
    properties
        rhat_C;  % Camera rays's directions.
        nr;  % Total number of rays (pixels). = size(1, rhat_C).
        p_BC_list; % = repmat(p_BC', 1, nr);  % As needed by OPCODE
        
        p_BC;
        tree;
    end
    
    properties (Constant)
        nh = 224;  % pixels in horizontal direction.
        nv = 171;  % pixels in vertical direction.
        fov_h = 62 * pi / 180;  % Horizontal FOV, rads.
        fov_v = 45 * pi / 180;  % Vertical FOV, rads.
        
    end
    
    methods
    
        % Constructor.
        % p_BC: position of the camera frame C in the bubble frame B.
        function this = PicoFlexCamera(p_BC, p_BP, t)
            % Position of the picoflex camera frame C in the bubble frame B.
            this.p_BC = p_BC;
            
            this.rhat_C = generate_picoflex_rays();
            this.nr = size(this.rhat_C, 1);
            this.p_BC_list = repmat(p_BC', 1, this.nr);  % As needed by OPCODE

            % Generate AABB tree for the mesh.
            addpath('../opcodemesh/matlab'); % Make OPCODE lib available.
            this.tree = opcodemesh(p_BP', t'); % NOTE!: I am using the transpose!
        end
        
        function [does_hit, dist, tri_index, bar_coos, p_BY] = GeneratePointCloud(this, p_BP, sigma_dist)
            
            this.tree.update(p_BP');
            [does_hit, dist, tri_index, bar_coos, p_BY] = ...
                this.tree.intersect(this.p_BC_list, this.rhat_C');
            p_BY = p_BY';

            % Add noise.
            %dist_noise = normrnd(0.0, sigma_dist, nr, 1); % You need a toolbox for
            %this! agghhh!
            dd = zeros(this.nr, 1);
            for ir = 1:this.nr
                dd(ir) = generate_gaussian_samples(0.0, sigma_dist);
                p_BY(ir, :) = p_BY(ir, :) + dd(ir) * this.rhat_C(ir, :);
            end
                        
        end
    
    end % End of methods section.
    
end