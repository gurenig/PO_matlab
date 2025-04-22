function c = center_slice(v, k)
    n = length(v);
    if k > n
        error('Requested window size k exceeds vector length.');
    end

    center = floor((n + 1)/2);  % Center index for both odd/even
    half_k = floor(k / 2);

    % Shift window so it's centered and has length k
    if mod(k, 2) == 1
        start_idx = center - half_k;
    else
        start_idx = center - half_k + 1;
    end

    end_idx = start_idx + k - 1;

    % Clamp to valid bounds
    if start_idx < 1
        start_idx = 1;
        end_idx = k;
    elseif end_idx > n
        end_idx = n;
        start_idx = n - k + 1;
    end

    c = v(start_idx:end_idx);
end
