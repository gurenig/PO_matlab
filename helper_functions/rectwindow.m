%> @file rectwindow.m
%> @brief Generates a rectangular window with tapered edges.
%>
%> Constructs a symmetric 1D window of length n, where the center region has
%> full amplitude (1) and the left and right edges are tapered down to a height h
%> over a width defined by d.
%>
%> @param d Width of taper transition (0 < d < 0.5)
%> @param h Edge height (0 < h < 1)
%> @param n Total number of samples in the window
%>
%> @retval val Output rectangular window as a 1D array
function val = rectwindow(d,h,n)
    % 0 < d < 0.5
    % 0 < h < 1
    x = linspace(-0.5, 0.5, n);
    val = ones(size(x));
    val(x < -0.5 + d) = h;
    val(x > 0.5 - d) = h;
end