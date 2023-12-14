# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 11:37:25 2023

@author: Mattc
"""

import matplotlib.pyplot as plt

def plot_data(file_path):
    # Read the data from the text file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract data from columns
    column1 = []
    column2 = []
    column3 = []

    for line in lines:
        data = line.split()
        column1.append(float(data[0]))
        column2.append(float(data[1]))
        column3.append(float(data[2]))

    # Plotting as scatter graph
    plt.scatter(column1, column2, label='Column 2', marker='o')
    plt.scatter(column1, column3, label='Column 3', marker='x')

    # Add labels and legend
    plt.xlabel('Column 1')
    plt.ylabel('Values')
    plt.legend()

    # Show the plot
    plt.show()

# Example usage
file_path = 'xf.txt'  # Replace with the path to your text file
plot_data(file_path)
