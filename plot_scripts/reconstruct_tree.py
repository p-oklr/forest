"""Reconstructs the tree representation from a tree description file."""

import graphviz
import numpy as np

DESC_FILE = "queue.txt"


def tree_recursive(dot, leaves, pbhs, current_number):
    """Explores the tree architecture and builds the graph.

    Parameters
    ==========
    - dot: object
        Graph object from graphviz
    - leaves: list of integer
        List of all the leaves in the description file
    - pbhs: list of integer
        List of all the pbhs in the description file
    - current_number: integer
        Integer ID that uniquely specifies a node
    """

    if not current_number in leaves:

        # You have two kids
        dot.node(f"{2 * current_number}")
        dot.node(f"{2 * current_number + 1}")

        dot.edge(f"{current_number}", f"{2 * current_number}")
        dot.edge(f"{current_number}", f"{2 * current_number + 1}")

        tree_recursive(dot, leaves, pbhs, 2 * current_number)
        tree_recursive(dot, leaves, pbhs, 2 * current_number + 1)

    # Node is dark if inside a PBH
    if current_number in pbhs:
        dot.node(f"{current_number}", style="filled", color="black", fontcolor="white")


def main(filename):
    """Main block of the code.

    Parameters
    ==========
    - filename: String
        Name of the file describing the tree
    """

    # Make the graph
    dot = graphviz.Digraph(filename + "_tree")

    # Read queue
    data = np.loadtxt(filename)

    # Build sets of leaves and PBHs
    leaves = [x for x, y in data if y != -1]
    pbhs = [x for x, y in data if y == -1]

    # Explore the tree and fill in the graph
    tree_recursive(dot, leaves, pbhs, 1)

    # Render the graph
    dot.render()


if __name__ == "__main__":
    main(DESC_FILE)
