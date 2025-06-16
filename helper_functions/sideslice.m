%> @file sideslice.m
%> @brief Extracts elements from both ends of a vector.
%>
%> Returns a new vector composed of the first floor(k/2) elements and the
%> last (k - floor(k/2)) elements from the input vector v. Useful for symmetric
%> slicing from the edges.
%>
%> @param v Input vector
%> @param k Total number of elements to extract (must be <= length(v))
%>
%> @retval s Vector containing the concatenated slices from both ends
function s = sideslice(v, k)
    n = length(v);
    if k > n
        error('k exceeds the length of the vector.');
    end

    left = floor(k / 2);
    right = k - left;

    s = [v(1:left), v(end-right+1:end)];
end
