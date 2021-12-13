import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

def get_job_id():
    return os.environ['PBS_JOBID'].split('.')[0]


def plot(name, df, x_axis, y_axis, category, plot_type='line', scale='linear', annotate=None):
  for idx, label in enumerate(df[category].unique()):
    sub_df = df[df[category] == label]
    color = plt.rcParams['axes.prop_cycle'].by_key()['color'][idx]
    if plot_type == 'scatter':
      plt.scatter(
        x=x_axis,
        y=y_axis,
        label=label,
        color=color,
        data=sub_df
      )
    else:
      plt.plot(
        x_axis,
        y_axis,
        label=label,
        color=color,
        data=sub_df
      )

    if annotate:
      for i in sub_df[annotate].index:
        plt.annotate(sub_df[annotate][i], (sub_df[x_axis][i], sub_df[y_axis][i]), color=color)

  plt.xlabel(x_axis)
  plt.ylabel(y_axis)
  plt.xscale(scale)
  plt.yscale(scale)
  plt.legend(title=category, loc='right')

  file_name = f'{name}.png'
  plt.savefig(file_name)
  plt.clf()
  print(f'Figure saved to {file_name}')

def plot_figure_one():
    name = f'figure_one_{get_job_id()}'
    df = pd.read_json(f'{name}.json')
    x_axis = 'n_asteroids'
    y_axis = 'asteroid time steps per second'
    category = 'precision'
    plot(name, df, x_axis, y_axis, category, scale='log', plot_type='scatter')

def plot_figure_two():
    name = f'figure_two_{get_job_id()}'
    df = pd.read_json(f'{name}.json')
    x_axis = 'simulated years per second'
    y_axis = 'error'
    category = 'precision'
    plot(name, df, x_axis, y_axis, category, scale='log', annotate='n_steps')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('figure')
    args = parser.parse_args()

    if args.figure == '1':
        plot_figure_one()
    elif args.figure == '2':
        plot_figure_two()




