import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})

def interpret_text_file(file_path):
    try:
        data = np.loadtxt(file_path, dtype=np.float64)

        col1 = data[:, 0]
        col2 = data[:, 1]
        col3 = data[:, 2]
        col4 = data[:, 3]

    except Exception as e:
        print(f"An error occurred: {e}")
        
    return col1, col2, col3, col4

file_path = "coords.txt"

c1, c2, c3, c4 = interpret_text_file(file_path)

x = np.concatenate((c3, -c3))
y = np.concatenate((c4,c4))


plt.scatter(x,y)
plt.show()