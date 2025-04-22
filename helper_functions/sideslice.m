function s = sideslice(v, k)
    n = length(v);
    if k > n
        error('k exceeds the length of the vector.');
    end

    left = floor(k / 2);
    right = k - left;

    s = [v(1:left), v(end-right+1:end)];
end
