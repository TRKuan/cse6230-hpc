import os

import matplotlib.pyplot as plt
import pandas as pd

def plot_figure():
    name = 'ring_test_169650'
    df = pd.read_json(f'{name}.json')
    x_axis = 'size'
    y_axis = 'seconds per message'

    plt.plot(
        x_axis,
        y_axis,
        data=df
    )
    plt.xlabel('number of processes')
    plt.ylabel(y_axis)

    file_name = f'{name}.png'
    plt.savefig(file_name)
    plt.clf()
    print(f'Figure saved to {file_name}')

if __name__ == '__main__':
    plot_figure()
