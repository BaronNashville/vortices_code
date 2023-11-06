function c = error_code_to_colour(err)
%Used for graph coloring

if err == 0
    c = [0 1 0];    %green      - Lyapunov stable
elseif err == 1
    c = [1 0 0];    %red        - unstable
elseif err == 2
    c = [1 0 1];    %magenta    - inconclusive test
elseif err == 3
    c = [1 1 0];    %yellow     - numercial uncertainty
elseif err == 4
    c = [0 0 0];    %black      - could not prove existence
elseif err == 5
    c = [0 0 1];    %blue       - used for testing
end
end