import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from numpy import array


def plot_height(p):
    CA = array([[.1, 51, 51, 153],
                [3.0, 0, 197, 120],
                [3.5, 0, 100, 120],
                [3.0, 0, 145, 120]])
    colors = CA[CA[:, 0].argsort()][:, 1:]/255.
    cmap = matplotlib.colors.ListedColormap(colors)
    bounds = [.1, 2.0, (int(p.hmap._data.max()) - 2) / 2 + 2,
              ((int(p.hmap._data.max()) - 2) * 3 / 4) + 2, int(p.hmap._data.max())]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    plt.imshow(p.hmap._data, extent=(0, p.wd.width - 1, p.wd.height - 1, 0),
               interpolation='nearest', cmap=cmap, norm=norm)
    plt.imshow(p.hmap._data, extent=(0, p.wd.width - 1, p.wd.height - 1, 0),
               interpolation='nearest', cmap=cmap, norm=norm)
    plt.show()
