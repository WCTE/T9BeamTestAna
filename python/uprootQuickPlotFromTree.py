import uproot
import matplotlib.pyplot as plt
import numpy as np

# jk 15.12.2023
# using chat gpt4

# Function to read a ROOT file, extract data, and create a 2D plot with a selection
def plot_2d_with_selection(root_file_path, tree_name, x_var, y_var, selection_var, selection_cut):
    # Open the ROOT file
    file = uproot.open(root_file_path)

    # Access the TTree
    tree = file[tree_name]

    # Extract data from the TTree as arrays of floats
    x_data = tree.arrays(x_var)
    y_data = tree.arrays(y_var)
    selection_data = tree.arrays(selection_var)

    # Apply selection cut
    #selected_indices = np.where(selection_data < selection_cut)
    #x_selected = x_data[:, selected_indices]
    #y_selected = y_data[:, selected_indices]

    x_selected = x_data
    y_selected = y_data

    
    # Create a 2D plot
    plt.scatter(x_selected[0], y_selected[0], s=5, alpha=0.5)
    plt.xlabel(f'{x_var}[0]')
    plt.ylabel(f'{x_var}[1]')
    plt.title(f'2D Plot with {selection_var} < {selection_cut}')
    plt.show()


# Example usage
root_file_path = "output/ntuple_000409.root"
tree_name = "ACT2L"
x_variable = "PeakVoltage"
y_variable = "PeakTime"
selection_variable = "nPeaks"
selection_cut_value = 2

plot_2d_with_selection(root_file_path, tree_name, x_variable, y_variable, selection_variable, selection_cut_value)
