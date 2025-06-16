%> @file coswindow.m
%> @brief Generates a custom cosine window with flattened sides.
%>
%> Creates a 1D symmetric window of length n that transitions using a cosine
%> shape over a width determined by d, and flattens to a base level of h at the ends.
%>
%> @param d Width of cosine transition (0 < d < 0.5)
%> @param h Floor level at edges (0 < h < 1)
%> @param n Total number of samples
%>
%> @retval val 1D array of length n containing the window values
function val = coswindow(d,h,n)
    % 0 < d < 0.5
    % 0 < h < 1
    x = linspace(-0.5, 0.5, n);
    val = (1 + cos((2 * pi * x) / (2 * d - 1))) * 0.5 * (1 - h) + h;
    val(x < -0.5 + d) = h;
    val(x > 0.5 - d) = h;
end
