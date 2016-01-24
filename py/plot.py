import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def plot_dists(dist, var_type, path):
    """Plot size density plot
    """
    fig, ax = plt.subplots()
    try:
        ax = sns.distplot(dist, hist=True)
        plt.title('{var_type} size distribution'.format(var_type=var_type))
        fig.savefig((path + '{var_type}_size_distribution.png'.
            format(var_type=var_type)))
    except:
        pass


