function p = bisection(f, a, b)
    if f(a)*f(b)>0 
        disp('Wrong choice bro')
        p = [];
    else
        p = (a + b)/2;
        err = abs(f(p));
        while err > 1e-6
            if f(a)*f(p)<0 
            	b = p;
            else
            	a = p;          
            end
            p = (a + b)/2; 
            err = abs(f(p));
        end
    end
end