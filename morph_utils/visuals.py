import numpy as np

def basic_morph_plot(morph,
                         ax,
                         title="",
                         morph_colors={1: 'black', 2: "firebrick", 4: "orange", 3: "steelblue"},
                         side=False,
                         xoffset=0,
                         line_w=0.5,
                         scatter=False,
                         scatter_roots=True,
                         dotsize=None,
                         scatter_soma = True,
                         soma_dot_size = 1,
                         plot_soma = True
                         ):
    """
    Populate a matplotlib axes with a morphology plot.

    :param morph: neuron_morphology Morphology object
    :param ax: matplotlib axes
    :param title: plot title
    :param morph_colors: dictionary to represent colors for each compartment
    :param side: if True will plot zy
    :param xoffset: How far to push the nodes in x. Useful when plotting multiple cells together
    :param line_w: linewidth for morphology lines
    :param scatter: if you want to scatter nodes in addition to drawing their lines
    :param scatter_roots: scatter nodes whose parent = -1
    :param dotsize: dot size for scatter points
    :param plot_soma: if True will plot lines connecting soma to its children. 
    :return: None, your axis will just be filled with neuron
    """
    ax.set_title(title)
    if dotsize is None:
        dotsize = line_w * 20
    if scatter:

        for ntype, color in morph_colors.items():
            if not side:
                nodes = np.array([[n['x'], n['y']] for n in morph.nodes() if n['type'] == ntype])
            else:
                nodes = np.array([[n['z'], n['y']] for n in morph.nodes() if n['type'] == ntype])

            try:
                ax.scatter(nodes[:, 0], nodes[:, 1], c=color, s=dotsize, label="Type {}".format(ntype))
            except:
                its = "missing a compartment"

    for compartment, color in morph_colors.items():
        lines_x = []
        lines_y = []
        for c in [n for n in morph.nodes() if n['type'] == compartment]:
            if c["parent"] == -1:
                continue
            try:
                p = morph.node_by_id(c["parent"])
            except:
                p = c

            if side:
                lines_x += [p["z"] + xoffset, c["z"] + xoffset, None]
            else:
                lines_x += [p["x"] + xoffset, c["x"] + xoffset, None]
            lines_y += [p["y"], c["y"], None]
        lines_y = [None if v is None else 1 * v for v in lines_y]
        ax.plot(lines_x, lines_y, c=color, linewidth=line_w, zorder=compartment)


    #scatter soma as per morph.get_soma()
    soma_node = morph.get_soma()
    if (soma_node) and (scatter_soma):
        if not side:
            ax.scatter(soma_node['x'] + xoffset, soma_node['y'], c='k', marker='X', s=soma_dot_size, label='morph.get_soma()')
        else:
            ax.scatter(soma_node['z'] + xoffset, soma_node['y'], c='k', marker='X', s=soma_dot_size, label='morph.get_soma()')

    #scatter root nodes throughout 
    if scatter_roots:
        root_nodes = [n for n in morph.nodes() if (n['parent'] == -1) and (n['type']==1)]
        for rn in root_nodes:
            if not side:
                ax.scatter(rn['x'] + xoffset, rn['y'], c='k', marker='X', s=dotsize * 2, label='root (parent=-1,type=1)')
            else:
                ax.scatter(rn['z'] + xoffset, rn['y'], c='k', marker='X', s=dotsize * 2, label='root (parent=-1,type=1)')

    soma_root = morph.get_soma()
    if (soma_root) and (plot_soma):
        for ch in morph.get_children(soma_root):
            if not side:
                ax.plot([ch['x'] + xoffset, soma_root['x']+ xoffset],
                        [ch['y'], soma_root['y']], c='k', linewidth=line_w * 2.5)
            else:
                ax.plot([ch['z'] + xoffset, soma_root['z']+ xoffset],
                        [ch['y'], soma_root['y']], c='k', linewidth=line_w * 2.5)
