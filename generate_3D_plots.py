import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Enable interactive mode for rotatable plots
plt.ion()  # Enables interactive mode

# Color mapping
color_mp = {
    'D': 'b', 'E': 'b',   # Blue for negative charge
    'R': 'r', 'K': 'r', 'H': 'r',  # Red for positive charge
    'N': 'y', 'Q': 'y', 'S': 'y', 'T': 'y', 'Y': 'y',  # Yellow for polar uncharged amino acids
    'A': 'g', 'V': 'g', 'L': 'g', 'I': 'g', 'P': 'g', 'F': 'g', 'M': 'g', 'W': 'g', 'C': 'g', 'G': 'g'  # Green for nonpolar amino acids
}

# Function to plot the 3D data
def plot_3d_data(data, title):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(title)

    for amino_acid, coordinates in data.items():
        x, y, z = coordinates
        color = color_mp.get(amino_acid, 'k')  # Default to black if not found in color_mp
        ax.scatter(x, y, z, c=color, label=amino_acid)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Use the legend with unique labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1.05, 1), ncol=2)

    # Add color key for the meanings of the colors
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Negative charge (D, E)', markerfacecolor='b', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Positive charge (R, K, H)', markerfacecolor='r', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Polar uncharged (N, Q, S, T, Y)', markerfacecolor='y', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Nonpolar (A, V, L, I, P, F, M, W, C, G)', markerfacecolor='g', markersize=10),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 0.7), ncol=1)

    plt.show()

# Your data
bindingsite_data = {
    'D': [-2.25, -0.7, -2.2],
    'E': [-2.4, -2.75, -2.4],
    'S': [-2.4, -2.75, -2.15],
    'G': [-2.15, -1.9, -2.7],
    'H': [-2.2, -2.0, -2.3],
    'F': [-2.2, -2.5, -2.0],
    'I': [-2.3, -2.7, -1.9],
    'K': [-2.25, -2.35, -2.6],
    'L': [-2.3, -3.2, -3.7],
    'M': [-2.35, -1.85, -3.4],
    'C': [-2.2, -2.1, -1.85],
    'P': [-2.7, -3.5, -3.3],
    'A': [-2.75, -3.15, -3.3],
    'N': [-2.25, -2.2, -1.9],
    'Q': [-2.3, -1.5, -1.9],
    'R': [-2.7, -3.3, -2.7],
    'T': [-2.8, -3.8, -2.0],
    'V': [-2.5, -3.75, -3.35],
    'W': [-2.7, -2.3, -1.95],
    'Y': [-3.2, -4.25, -2.4],
    'Z': [-2.8, -3.1, -1.5],
    'U': [-4.0, -3.0, -1.5],
    'B': [-3.75, -3.45, -1.6]
}

alternativesplicing_data = {
    'W': [-3.45, -1.25, -2.7],
    'H': [-3.5, -1.45, -2.3],
    'Y': [-3.25, -1.65, -2.8],
    'Q': [-2.95, -1.65, -2.85],
    'G': [-2.6, -1.7, -2.45],
    'I': [-3.25, -1.8, -2.5],
    'C': [-3.0, -1.8, -2.45],
    'R': [-2.45, -1.95, -2.3],
    'S': [-2.55, -2.05, -2.9],
    'L': [-2.75, -2.2, -2.4],
    'V': [-2.65, -2.15, -2.4],
    'D': [-2.45, -2.2, 0.5],
    'F': [-2.25, -2.25, -2.25],
    'M': [-2.05, -2.4, -2.15],
    'E': [-1.7, -2.3, 0.5],
    'T': [-1.55, -2.45, -1.95],
    'A': [-1.8, -2.6, 2.2],
    'P': [-1.75, -2.65, 2.25],
    'K': [-1.5, -2.6, -2.25],
    'U': [-3.5, -2.2, -1.6],
    'N': [0.0, -2.75, -2.0]
}

threeDprotein_data = {
    'L': [-2.75, -0.2, -2.7],
    'R': [-0.9, -1.25, -1.0],
    'K': [-0.95, -1.25, -0.6],
    'P': [-1.0, -1.2, -0.75],
    'V': [-1.1, -1.2, -0.8],
    'I': [-1.1, -1.2, -0.85],
    'D': [-1.0, -1.4, -2.5],
    'N': [-1.0, -1.45, 0.0],
    'Y': [-0.95, -1.4, -0.0],
    'U': [-1.0, -1.6, -2.5],
    'A': [-1.0, -1.25, -0.75],
    'S': [-1.2, -1.4, -0.0],
    'T': [-1.3, -1.45, 0.05],
    'F': [-1.45, -1.4, -0.8],
    'M': [-1.5, -1.3, -0.85],
    'G': [-1.25, -1.6, 0.8],
    'Q': [-1.5, -1.45, 0.0],
    'C': [-1.8, -1.45, -0.75],
    'W': [-2.05, -1.4, -0.85],
    'E': [-2.05, -1.6, -2.5],
    'H': [-2.7, -1.6, -0.55],
    'Z': [-2.85, -1.4, -2.45]
}

# Plot all three datasets
plot_3d_data(bindingsite_data, 'Binding Site Data')
plot_3d_data(alternativesplicing_data, 'Alternative Splicing Data')
plot_3d_data(threeDprotein_data, '3D Protein Data')

# Keep the plots interactive
plt.ioff()  # Turn off interactive mode if needed later
input("Press [Enter] to close the plots and exit...")
