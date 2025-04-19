function cmap = redbluecmap(m)
    % redbluecmap  - Diverging colormap from blue to white to red
    % Usage:
    %   colormap(redbluecmap(256))

    if nargin < 1
        m = 256;  % default colormap size
    end

    bottom = [0 0 1];     % Blue
    middle = [1 1 1];     % White
    top = [1 0 0];        % Red

    % Number of colors in each half
    if mod(m, 2) == 0
        n = m / 2;
        r = [linspace(bottom(1), middle(1), n), linspace(middle(1), top(1), n)];
        g = [linspace(bottom(2), middle(2), n), linspace(middle(2), top(2), n)];
        b = [linspace(bottom(3), middle(3), n), linspace(middle(3), top(3), n)];
    else
        n = floor(m / 2);
        r = [linspace(bottom(1), middle(1), n), middle(1), linspace(middle(1), top(1), n)];
        g = [linspace(bottom(2), middle(2), n), middle(2), linspace(middle(2), top(2), n)];
        b = [linspace(bottom(3), middle(3), n), middle(3), linspace(middle(3), top(3), n)];
    end

    cmap = [r(:), g(:), b(:)];
end
