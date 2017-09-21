function [fx,xr] = NewtRaph(x0,err,a,b,c)
    i = 0;
    xold = x0;
    x = x0;
    while(1)
        i = i + 1;
        fx = lab1fn_Imp(x,a,b,c);
        dfx = funct_deriv(x,a,b);
        x = x - fx/dfx;
        ea = abs((x-xold)/x)*100;
        if ea < err
            xr = x;
            break
        end
        xold = x;
    end
end
        