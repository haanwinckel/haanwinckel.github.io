using Random, Dates, Statistics, Plots, FileIO, Images

# Set a random seed based on the current time.
Random.seed!(Dates.value(Dates.now()))

#Will do everything twice, one for header and one for favicon 
for favicon in [false, true]

    if favicon
        width = 32      # output width in pixels
        height = 32      # output height in pixels
    else
        # High-resolution parameters
        width = 3360      # output width in pixels
        height = 800      # output height in pixels
    end

        

    # -------------------------------------------------------------------
    # Functions for finding an interesting Mandelbrot window
    # -------------------------------------------------------------------

    # Common Mandelbrot function
    function mandelbrot(c, max_iter)
        z = 0 + 0im
        n = 0
        while abs(z) ≤ 2 && n < max_iter
            z = z*z + c
            n += 1
        end
        return n
    end

    # Compute a low-resolution Mandelbrot grid for a given window.
    function compute_mandelbrot_window(x_min, x_max, y_min, y_max, width, height, max_iter)
        xs = range(x_min, x_max, length=width)
        ys = range(y_min, y_max, length=height)
        grid = [mandelbrot(xs[j] + ys[i]*im, max_iter) for i in 1:height, j in 1:width]
        return grid
    end

    # Compute the Shannon entropy of the grid.
    function measure_entropy(grid, max_iter)
        data = vec(grid)
        counts = zeros(Float64, max_iter + 1)
        for n in data
            counts[n+1] += 1
        end
        total = sum(counts)
        p = counts / total
        p_nonzero = filter(x -> x > 0, p)
        return -sum(p_nonzero .* log.(p_nonzero))
    end

    # Compute an average gradient (difference) measure.
    function measure_gradient(grid)
        dx = abs.(diff(grid, dims=2))
        dy = abs.(diff(grid, dims=1))
        return mean(dx) + mean(dy)
    end

    # Decide whether the candidate window is interesting.
    function is_interesting(grid, max_iter; 
            entropy_thresh=favicon ? 5.0 : 3.0, 
            gradient_thresh=favicon ? 0.75 : 0.5, 
            low_thresh=favicon ? 0.25 : 0.1, 
            high_thresh=favicon ? 0.75 : 0.9)
        count_diverge = count(n -> n < max_iter, grid)
        total = length(grid)
        frac = count_diverge / total
        if low_thresh < frac < high_thresh
            ent = measure_entropy(grid, max_iter)
            grad = measure_gradient(grid)
            # Accept if either the entropy or gradient measure is high enough.
            return ((ent > entropy_thresh) || (grad > gradient_thresh), max(ent/entropy_thresh, grad/gradient_thresh))
        else
            return (false, -Inf)
        end
    end

    # Loop to search for an interesting window.
    function find_interesting_window(max_attempts=100000)
        max_iter = 300
        best_distance = -Inf
        best_window = (Inf, Inf, Inf, Inf)
        for attempt in 1:max_attempts
            # Define ranges for a candidate center (you can adjust these ranges)
            #x_low, x_high, y_low, y_high = -1.34, -1.30, -0.12, -0.06
            x_low, x_high, y_low, y_high = -2.0, 0.5, 0.0, 1.0
            center = x_low + rand()*(x_high - x_low) + im*(y_low + rand()*(y_high - y_low))
            # Randomly choose a zoom factor (smaller means more zoomed in)
            zoom = 10.0^rand(-10:-8)
            x_offset = zoom * (x_high - x_low)
            
            x_min = real(center) - x_offset
            x_max = real(center) + x_offset
            # Maintain an aspect ratio corresponding to 3360 x 800 (i.e. height/width = 800/3360)
            ratio = height / width
            y_min = imag(center) - x_offset * ratio
            y_max = imag(center) + x_offset * ratio
            
            # Compute a low-res grid (for speed, e.g., 200x67 pixels)
            grid = compute_mandelbrot_window(x_min, x_max, y_min, y_max, min(200, width), min(67, height), max_iter)
            (good, value) = is_interesting(grid, max_iter)
            if good
                println("Found interesting window on attempt $attempt. Zoom was $zoom.")
                return (x_min, x_max, y_min, y_max)
            else
                if value > best_distance
                    best_distance = value
                    best_window = (x_min, x_max, y_min, y_max)
                end
            end
        end
        println("Failed to find an interesting window within $max_attempts attempts; using best candidate.")
        return best_window
    end

    # -------------------------------------------------------------------
    # Use the above functions to find a window.
    # -------------------------------------------------------------------
    window = find_interesting_window()
    println("Selected window: ", window)
    (x_min, x_max, y_min, y_max) = window

    # -------------------------------------------------------------------
    # High-resolution fractal image generation
    # -------------------------------------------------------------------

    max_iter = 300

    # Create high-res coordinate arrays based on the found window.
    xs = range(x_min, x_max, length=width)
    ys = range(y_min, y_max, length=height)
    mandelbrot_set = zeros(Int, height, width)
    for j in 1:width
        for i in 1:height
            c = xs[j] + ys[i]*im
            mandelbrot_set[i, j] = mandelbrot(c, max_iter)
        end
    end

    # Define a colormap with mostly white background and subtle light-blue fractal details.
    if favicon
        colors = ["#ffffff", "#b0c4de", "#4682b4", "#132265"]
    else
        colors = ["#ffffff", "#f0faff", "#e0f5ff", "#cceeff"]
    end
    cmap = cgrad(colors, max_iter)

    border_size = favicon ? 13 : 60

    # Generate the heatmap using bilinear interpolation and high DPI.
    heatmap(xs, ys, mandelbrot_set,
        size = (width+2*border_size, height+2*border_size),  # Force the output size to 3360 x 800 pixels
        xlims = (x_min, x_max),
        ylims = (y_min, y_max),
        color = cmap,
        axis = nothing,
        legend = false,
        interpolation = :bilinear, #dpi = 300,
        widen = false
    )

    filename = favicon ? "favicon.png" : "header_background.png"
    savefig(filename)

    # Optionally load the image and crop extra borders.
    img = load(filename)
    h, w = size(img)
    println("Generated image dimensions: ", (h, w))
    # Here we crop 60 pixels from each side (adjust as needed)
    cropped_img = img[(1+border_size):(h-border_size+1), (1+border_size):(w-border_size+1)]
    save(filename, cropped_img)

end