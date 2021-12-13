import matplotlib.pyplot as plt
import pandas as pd


if __name__ == '__main__':
    file_name = 'buffer.json'
    file_var_name = 'vector.json'

    df = pd.read_json(file_name)
    df_var = pd.read_json(file_var_name)

    x_axis = 'number of grid points per direction'
    y_axis = 'Mean lattice updates per second'

    plt.plot(
        x_axis,
        y_axis,
        label="packing and unpacking send and receive buffers",
        data=df
    )

    plt.plot(
        x_axis,
        y_axis,
        label="in place buffers with MPI vector types",
        data=df_var
    )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.legend()

    plt.savefig('figures.png')