using Plots

# Function to generate random x and y values for demonstration purposes
generate_random_data(n) = (rand(n), rand(n))

# Number of plots to generate
num_plots = 7

# Initialize an empty array to hold the individual plots
individual_plots = []

# Generate and add the plots to the array
for i in 1:num_plots
    x, y = generate_random_data(10)  # Generate 10 random data points for each plot
    p = Plots.plot(x, y, title = "Plot $i", xlabel = "x", ylabel = "y", legend = false)
    push!(individual_plots, p)
end

# Combine all individual plots into a single grid layout
Plots.plot(individual_plots..., layout = (2, 4), size = (1200, 600))
