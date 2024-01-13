function [controlInput,flags] = nmpcControl(i,x1,x2,controlInput,K,P,alpha)
    if [x1(i);x2(i)]'*P*[x1(i);x2(i)] > alpha
        [u, flags(i)] = ocp_van_der_pol(x1(i), x2(i), P, alpha);
        if isnan(u)
            if i == 1
                % catch if initial problem is infeasible
                controlInput(i) = 0;
            else
                % Handling of infeasible solutions
                controlInput(i) = controlInput(i-1);
                % controlInput(i) = -controlInput(i-1);
                % controlInput(i) = 0;
                % controlInput(i) = -K*[x1(i);x2(i)];
            end
        else
            % normal MPC iteration
            controlInput(i) = u;
        end
    else
        controlInput(i) = -K*[x1(i);x2(i)];
        flags(i) = 0;
    end
end