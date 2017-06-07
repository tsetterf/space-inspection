classdef InertiaRatios < handle
    %INERTIARATIOS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gamma = [1.1 1.05]';          % Inertia ratios gamma1 = J1/J3, gamma2 = J2/J3
        y = [0.05 0.05]';             % Equivalent values of gamma in intermediate space Y
    end
    
    methods
        %% Constructor
        function this = InertiaRatios(gamma)           
            this.setGamma(gamma);
        end
        %% Getters and setters
        function g = getGamma(this)
            g = this.gamma;
        end
        function setGamma(this, g)    
            this.gamma = g;
            % Update equivalent of gamma in intermediate space Y
            this.y = [g(1) - g(2); g(2) - 1];
        end
        %% Retract from R^2 to space of valid inertia ratios Gamma
        function ratios = retract(this, x)
            % Get intermediate mapping from R^2 to Y
            y  = [this.logisticFunction(x(1)); exp( x(2)/this.y(2) + log(this.y(2)) ) ];

            % Get final mapping from Y to Gamma
            g = [y(1) + y(2) + 1; y(2) + 1];
            
            % Return inertia ratios object
            ratios = InertiaRatios(g);
        end
        %% Get local coords (R^2) of a variable by lifting g from space of valid inertia ratios to R^2
        function x = localCoordinates(this, ratios)
            
            % Get inertia ratios gamma
            g = ratios.getGamma();
            
            % Get the intermediate mapping from Gamma to Y
            y = [g(1) - g(2); g(2) - 1];

            % Get the final mapping from Y to R^2
            x = [this.inverseLogisticFunction(y(1)); this.y(2) * (log(y(2)) - log(this.y(2)))];
            
        end
        %% Get the matrix form of the inertia tensor
        function J = matrix(this)
            J = diag([this.gamma(1) this.gamma(2) 1]);
        end
        %% Get the inverse of the matrix form of the inertia tensor
        function Jinv = matrixInverse(this)
            Jinv = diag([1/this.gamma(1) 1/this.gamma(2) 1]);
        end
        %% Print the ratios
        function print(this)
            disp(['gamma: ' num2str(this.gamma(1)) ' ' num2str(this.gamma(2))]);
        end
        %% Check if the ratios are equal
        function isEqual = equals(this, ratios, tol)
            if (abs(this.gamma(1) - ratios.gamma(1)) < tol) && (abs(this.gamma(2) - ratios.gamma(2)) < tol)
                isEqual = 1;
            else
                isEqual = 0;
            end
        end
        %% The sigmoid function for mapping from x1 (in R^2) to y1 (in Y)
        function y1 = logisticFunction(this, x1)
            % Get the values of the generalized logistic function sigmoid over t given a y intercept of yb. 
            % See logbook #2 pp. 168 and https://en.wikipedia.org/wiki/Generalised_logistic_function, and 
            % physical formula sheet.
            % y(t) = A +          K - A
            %            --------------------------
            %            ( C + Q*exp(-B*x1) )^(1/nu)
            % Lots of sigmoids defined here too, in case this doesn't work out:
            % https://www.researchgate.net/profile/Even_Tjorve/publication/227734985_Shapes_and_functions_of_speciesarea_curves_a_review_of_possible_models/links/54ba83b60cf24e50e9403313.pdf
            
            % Choose basic parameters
            A = 0;
            K = 1;
            C = 1;
            yb = this.y(1);
            
            % Deal with unsolvable cases when yb/K < 0.5
            if yb/K <= 0.5
                yb = (1-yb/K) * K;
                useSymmetry = 1;
            else
                useSymmetry = 0;
            end

            % Solve 1 = (yb/K)^nu * (nu+1) to have d^2y/dt^2 = 0 at yb
            nu = lambertw(-1,yb/K * log(yb/K)) / log(yb/K) - 1;

            % Sub in nu to B = nu / (1 - (yb/K)^nu * yb) to have dy/dt = 1 at yb
            B =  nu / ( (1 - (yb/K)^nu) * yb);

            % Compute Q such that y(0) = yb
            Q = -1 + (K/yb)^nu;

            % Compute function
            if useSymmetry == 0
                y1 = A + (K-A) ./ ( (C + Q*exp(-B*x1)).^(1/nu) );
            else
                y1 = K - (K-A) ./ ( (C + Q*exp( B*x1)).^(1/nu) );
            end

        end
        %% The inverse of the sigmoid function. Maps y2 (in Y) to x2 (in R^2)
        function x1 = inverseLogisticFunction(this, y1)
            % Get the values of the inverse of the generalized logistic function sigmoid given values y and a 
            % y-intercept of yb. 
            % See logbook #2 pp. 168 and https://en.wikipedia.org/wiki/Generalised_logistic_function, and 
            % physical formula sheet.
            % Formulas below are for cases where yb > 1/exp(1)
            % y(t) = A +          K - A
            %            --------------------------
            %            ( C + Q*exp(-B*t) )^(1/nu)
            %                   / (K/y)^nu - 1 \
            % t(y) = -1/B * log| -------------- |
            %                   \      Q       /
            % Lots of sigmoids defined here too, in case this doesn't work out:
            % https://www.researchgate.net/profile/Even_Tjorve/publication/227734985_Shapes_and_functions_of_speciesarea_curves_a_review_of_possible_models/links/54ba83b60cf24e50e9403313.pdf

            % Choose basic parameters
            A = 0;
            K = 1;
            C = 1;
            yb = this.y(1);

            % Deal with unsolvable cases when yb/K < 1/exp(1)
            if yb/K <= 0.5
                yb = (1-yb/K) * K;
                useSymmetry = 1;
            else
                useSymmetry = 0;
            end

            % Solve 1 = (yb/K)^nu * (nu+1) to have d^2y/dt^2 = 0 at yb
            nu = lambertw(-1,yb/K * log(yb/K)) / log(yb/K) - 1;

            % Sub in nu to B = nu / (1 - (yb/K)^nu * yb) to have dy/dt = 1 at yb
            B =  nu / ( (1 - (yb/K)^nu) * yb);

            % Compute Q such that y(0) = yb
            Q = -1 + (K/yb)^nu;

            % Compute function
            if useSymmetry == 0
                x1 = -1/B * log( ((K./y1).^nu     - 1) / Q);
            else
                x1 =  1/B * log( ((K./(K-y1)).^nu - 1) / Q);
            end

        end
    end
    
end

