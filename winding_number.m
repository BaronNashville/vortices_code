function winding = winding_number(A, corners, num_steps)
%A: matrix we are trying to verify eigenvalues for
%corners: array containing the 4 corners of our contour
%num_steps: the number of steps taken along each segment

N = size(A,1);
M = size(A,2);

if N ~= M
    winding = NaN;
else
    %Used for less confusing code later on
    B = @(z) A - z*eye(N);
    
    %Eigenvalue equation
    g = @(z) det_intval(B(z));
    
    %Initializing variables
    m = length(corners);    
    err = 0;
    
    possible_crossings = zeros(m,num_steps);
    derivative_dir = zeros(m,num_steps);
    
    points = zeros(m,num_steps);
    radii = zeros(m,num_steps);
    
    for i = 1:m
        start = corners(i);
        finish = corners(1+mod(i,m));
        
        %Derivative of gamma wrt t
        dgamma = finish-start;
        
        %Step size
        delta_z = dgamma/(num_steps);
        
        %Segment evaluated at partition point s
        gamma = @(s) start + s*delta_z;

        for n = 1:num_steps
            %Ball containing the line segment between gamma(n-1) and gamma(n)
            x = hull(gamma(n-1),gamma(n));
            
            tmp = g(x);
            
            if isnan(tmp)
                err = 1;
                break
            end
            
            c = 'b';
            points(i,n) = mid(tmp); radii(i,n) = rad(tmp);
            
            if (~(imag(tmp) < 0)  && ~(imag(tmp) > 0)) && (~(real(tmp) < 0)  && ~(real(tmp) > 0))
                %The number could be 0 and then we run into issues
                err = 1;
            elseif (~(imag(tmp) < 0)  && ~(imag(tmp) > 0)) && real(tmp) > 0
                %We may have crossed the ray with angle 0
                c = 'r';
                possible_crossings(i,n) = 1;
            end
            
            plotintval(tmp, c);
            hold on
            
            %Imaginary component of the derivative of g along the current segment
            Im_dir = imag(g(x) * -dgamma * trace(B(x)^-1));
            
            if Im_dir > 0
                derivative_dir(i,n) = 1;
            elseif Im_dir < 0
                derivative_dir(i,n) = -1;
            end
            
        end
        
    end
    
    winding = 0;
    
    if err
        winding = NaN;
    else
        for i = 1:m
            for n = 1:num_steps
                if possible_crossings(i,n) == 1
                    if n == 1
                        %If this is the first step, we include it in the calculation
                        winding = winding + derivative_dir(i,n);
                    else
                        if derivative_dir(i,n) == 0
                            %Parrallel crossing
                            err = 1;                            
                        elseif ~(possible_crossings(i,n-1) == 1 && derivative_dir(i,n-1) == derivative_dir(i,n))
                            %We only include it in the calculation if we haven't already included this crossing
                            winding = winding + derivative_dir(i,n);
                        end
                    end
                end  
            end
        end
    end
    
    if err
        winding = NaN;
    end
    
    max_x = max(max(real(points)));
    plot(0:ceil(max_x),zeros(1,ceil(max_x)+1), 'k')
end