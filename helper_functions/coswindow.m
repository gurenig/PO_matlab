function val = coswindow(d,h,n)
    % 0 < d < 0.5
    % 0 < h < 1
    x = linspace(-0.5,0.5,n);
    val = (1+cos((2*pi*x)/(2*d-1)))*0.5*(1-h)+h;
    val(x<-0.5+d) = h;
    val(x>0.5-d) = h;
end

