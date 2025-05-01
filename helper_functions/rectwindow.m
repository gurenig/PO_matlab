function val = rectwindow(d,h,n)
    % 0 < d < 0.5
    % 0 < h < 1
    x = linspace(-0.5,0.5,n);
    val = ones(size(x));
    val(x<-0.5+d) = h;
    val(x>0.5-d) = h;
end

