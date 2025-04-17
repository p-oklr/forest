"""Reconstructs the 2D map representation from a tree description file."""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

plt.rc("text", usetex=True)
plt.rc("font", family="serif")

MAX_WIDTH = 1
MAX_HEIGHT = 1
DESC_FILE = "queue.txt"


class MapObject:
    """Utilities to make a map representation of a single tree"""

    def __init__(self, ax, data):
        """Initializes the map object

        Parameters
        ==========
        - ax: Matplotlib axis
            Output axis for the map
        - data: numpy array
            Data describing the tree
        """

        self.ax = ax
        self.data = data

        # Build sets of leaves and PBHs
        self.leaves = [x for x, y in data if y != -1]
        self.pbhs = [x for x, y in data if y == -1]

        # Find the range for the number of efolds
        self.min_efold, self.max_efold = min(data[:, 1]), max(data[:, 1])

    def get_color(self, efold):
        """Assigns a color in the colormap.

        Parameters
        ==========
        - efold: float
            Number of efolds

        """
        return colormaps["Oranges"](
            (efold - self.min_efold) / (self.max_efold - self.min_efold)
        )

    def make_map(
        self,
        current_number,
        x0,
        y0,
        width,
        height,
    ):
        """Builds a map .

        Parameters
        ==========
        - current_number: integer
            Integer ID that uniquely specifies a node
        - x0: float
            Leftmost edge of the box
        - y0: float
            Bottom edge of the box
        - width: float
            Width of the box
        - height: float
            Height of the box
        """

        if current_number in self.leaves:
            efold = self.data[np.where(self.data[:, 0] == current_number)[0], 1]

            self.ax.add_artist(
                mpatches.Rectangle(
                    (x0, y0),
                    width,
                    height,
                    facecolor=self.get_color(efold),
                )
            )

        else:

            # Split axis determined by
            if width >= height:
                # Split along x
                self.make_map(
                    2 * current_number,
                    x0,
                    y0,
                    width / 2,
                    height,
                )
                self.make_map(
                    2 * current_number + 1,
                    x0 + width / 2,
                    y0,
                    width / 2,
                    height,
                )

            else:
                # Split along y
                self.make_map(
                    2 * current_number,
                    x0,
                    y0,
                    width,
                    height / 2,
                )
                self.make_map(
                    2 * current_number + 1,
                    x0,
                    y0 + height / 2,
                    width,
                    height / 2,
                )

        if current_number in self.pbhs:
            self.ax.add_artist(
                mpatches.Rectangle(
                    (x0, y0), width, height, facecolor="k", edgecolor="white"
                )
            )


def main(filename):
    """Main block of the code.

    Parameters
    ==========
    - filename: String
        Name of the file describing the tree
    """

    # Read queue
    data = np.loadtxt(filename)

    ########################################
    #           Build the figure           #
    ########################################

    # Define the axis
    fig, ax = plt.subplots(figsize=(4, 3))

    # Instanciate the map
    tree_map = MapObject(ax, data)

    # Explore the tree and fill in the map
    tree_map.make_map(
        1,
        0,
        0,
        MAX_WIDTH,
        MAX_HEIGHT,
    )

    # create dummy invisible image
    # (use the colormap you want to have on the colorbar)
    img = plt.imshow(
        np.array([[tree_map.min_efold, tree_map.max_efold]]), cmap="Oranges"
    )
    img.set_visible(False)
    plt.colorbar(orientation="vertical", label=r"$\mathcal{N}_{1\to j}$")

    # Axis properties
    ax.set_xlim(0, MAX_WIDTH)
    ax.set_ylim(0, MAX_HEIGHT)
    ax.set_aspect("equal", "box")
    ax.set_axis_off()
    fig.tight_layout()

    # Save the figure
    fig.savefig(filename + "_pattern.pdf", transparent=True)


if __name__ == "__main__":
    main(DESC_FILE)
