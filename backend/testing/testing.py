import pandas as pd
import plotly.express as px

def plot_umap_and_obs(umap_coords, data_df, color_column=None, template="plotly_white"):
    """
    Plot UMAP scatter plot using Plotly based on colors from a specified column.

    Parameters:
    - umap_coords: numpy array containing UMAP coordinates.
    - data_df: DataFrame containing metadata and columns for coloring.
    - color_column: Column name in data_df to color the scatter plot.
    - template: Plotly template for styling the plot.
    """
    
    # Convert UMAP coordinates to DataFrame
    umap_df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
    
    # Concatenate UMAP DataFrame with data_df side by side
    merged_df = pd.concat([umap_df, data_df.reset_index(drop=True)], axis=1)
    num_unique_classes = merged_df[color_column].nunique()

    # Create the scatter plot
    fig = px.scatter(
        merged_df,
        x='UMAP1',
        y='UMAP2',
        color=color_column if num_unique_classes < 30 else None,
        color_continuous_scale=None if num_unique_classes < 30 else 'Viridis',
        template=template
)

    show_legend = True if num_unique_classes < 30 else False

    # Update the layout and settings for the plot
    fig.update_layout(
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            linecolor='gray',
            linewidth=1.5,
            mirror=True,
            showline=True,
            showspikes=False,  # Hide axis spikelines
            title_text=""  # Empty title to hide it
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            linecolor='gray',
            linewidth=1.5,
            mirror=True,
            showline=True,
            showspikes=False,  # Hide axis spikelines
            title_text=""  # Empty title to hide it
        ),
        showlegend=show_legend,
        margin=dict(l=40, r=40, b=60, t=60),  # Adjusted bottom margin
        plot_bgcolor='white',
        width=600,
        height=450,
        annotations=[
            # X-axis label
            dict(
                xref="paper",
                yref="paper",
                x=0.5,
                y=-0.1,  # Adjusted y position
                showarrow=False,
                text="UMAP1",
                font=dict(size=14)
            ),
            # Y-axis label
            dict(
                xref="paper",
                yref="paper",
                x=-0.08,
                y=0.5,
                showarrow=False,
                text="UMAP2",
                textangle=-90,
                font=dict(size=14)
            )
        ]
    )
    

    # Display the plot
    fig.update_traces(marker=dict(size=4))
    fig.show()

# Sample call to the function (commented out for execution here)
# plot_umap_scatter(adata.obsm["X_umap"], adata.obs, color_column='leiden')


def add_dendrogram_to_plot(fig, adata, color_mapping=None):
    """
    Add a dendrogram to an existing Plotly figure.
    
    This function extracts dendrogram information from the AnnData object 
    and plots it onto the provided figure. The dendrogram typically represents
    the hierarchical clustering structure of the data.

    Parameters:
    ----------
    fig: plotly.graph_objects.Figure
        The Plotly figure to which the dendrogram will be added.
        
    adata: AnnData
        An AnnData object, which is expected to contain dendrogram data 
        under the 'dendrogram_leiden' key in its uns attribute.
        
    color_mapping: dict, optional
        A dictionary mapping the original colors in the dendrogram to desired 
        colors for visualization. If a color is not provided in this mapping,
        it defaults to black.

    Returns:
    -------
    None
    """
    # Extract dendrogram information from the AnnData object
    icoord = adata.uns["dendrogram_leiden"]['dendrogram_info']['icoord']
    dcoord = adata.uns["dendrogram_leiden"]['dendrogram_info']['dcoord']
    color_list = adata.uns["dendrogram_leiden"]['dendrogram_info']['color_list']

    # Plot dendrogram lines onto the figure
    for xs, ys, color in zip(icoord, dcoord, color_list):
        line_color = color_mapping.get(color, 'black')
        fig.add_trace(
            go.Scatter(
                x=ys,
                y=xs,
                mode='lines',
                line=dict(color=line_color, width=1),
                hoverinfo='none'
            ),
            row=1, col=2
        )

        
def prepare_data(adata, cell_types, n_genes, gene_list, group_by, flavor, min_logfc=1):
    """
    Prepare gene expression data for plotting. 

    This function processes the AnnData object to extract the mean expression 
    and fraction of cells expressing specific genes, either provided through 
    a gene list or by identifying top expressed genes.

    Parameters:
    ----------
    adata : AnnData
        An object containing gene expression data and potentially 
        differential expression results.
        
    cell_types : list
        List of cell types to consider.
        
    n_genes : int
        Number of top genes to consider if gene_list is not provided.
        
    gene_list : list, optional
        List of specific genes to consider. If provided, n_genes is ignored.
        
    group_by : str
        Column name in adata.obs that contains cell type information.
        
    flavor : str
        Specifies the method for obtaining gene values. Should be one of 
        'mean_expression' or 'fold_change'. If 'mean_expression', the function 
        will retrieve the mean expression of each gene. If 'fold_change', 
        the function will fetch the fold change value from differential 
        expression results stored in adata.uns.
        
    min_logfc : float, optional
        Minimum log fold change to consider when selecting genes based 
        on fold change. Only used if flavor is 'fold_change'.

    Returns:
    -------
    tuple
        A tuple containing three arrays:
        - Mean expression data across cell types for selected genes.
        - Fraction of cells expressing each gene across cell types.
        - List of all genes considered.
    """
    data = []
    fraction_expressing_data = []
    all_genes = []
    
    if gene_list is not None:  # If gene_list is provided, collect data once for each gene in gene_list
        selected_genes = gene_list
        all_genes.extend(selected_genes)
        
        for gene in selected_genes:
            values = []
            fraction_values = []
            
            for ct in cell_types:
                subset = adata[adata.obs[group_by] == ct]
                
                if flavor == "mean_expression":
                    value = subset[:, gene].X.mean()
                elif flavor == "fold_change":
                    value = adata.uns["rank_genes_groups"]["logfoldchanges"][ct][
                        np.where(adata.uns["rank_genes_groups"]["names"][ct] == gene)
                    ][0]
                else:
                    raise ValueError("Invalid flavor. Choose 'mean_expression' or 'fold_change'.")
                
                fraction_value = np.sum(subset[:, gene].X > 0) / subset.shape[0]
                
                values.append(value)
                fraction_values.append(fraction_value)
            
            data.append(values)
            fraction_expressing_data.append(fraction_values)
    
    else:  # Original logic when gene_list is not provided
        for cell_type in cell_types:
            top_genes = adata.uns["rank_genes_groups"]["names"][cell_type][:n_genes * 2].tolist()  # Taking 2x to have a buffer
            selected_genes = []
            
            for gene in top_genes:
                if len(selected_genes) >= n_genes:
                    break  # Stop if we already have enough genes
                
                if flavor == "fold_change":
                    logfc_value = adata.uns["rank_genes_groups"]["logfoldchanges"][cell_type][
                        np.where(adata.uns["rank_genes_groups"]["names"][cell_type] == gene)
                    ][0]
                    
                    if logfc_value >= min_logfc:
                        selected_genes.append(gene)
                else:
                    selected_genes.append(gene)
            
            # Check if we have enough genes; if not, fill with top genes irrespective of fold change
            remaining_genes_needed = n_genes - len(selected_genes)
            for gene in top_genes:
                if gene not in selected_genes:
                    selected_genes.append(gene)
                    remaining_genes_needed -= 1
                
                if remaining_genes_needed <= 0:
                    break
            
            all_genes.extend(selected_genes)
            
            for gene in selected_genes:
                values = []
                fraction_values = []
                
                for ct in cell_types:
                    subset = adata[adata.obs[group_by] == ct]
                    
                    if flavor == "mean_expression":
                        value = subset[:, gene].X.mean()
                    elif flavor == "fold_change":
                        value = adata.uns["rank_genes_groups"]["logfoldchanges"][ct][
                            np.where(adata.uns["rank_genes_groups"]["names"][ct] == gene)
                        ][0]
                    else:
                        raise ValueError("Invalid flavor. Choose 'mean_expression' or 'fold_change'.")
                    
                    fraction_value = np.sum(subset[:, gene].X > 0) / subset.shape[0]
                    
                    values.append(value)
                    fraction_values.append(fraction_value)
                
                data.append(values)
                fraction_expressing_data.append(fraction_values)
    
    return np.array(data), np.array(fraction_expressing_data), all_genes


def add_annotations(fig, cell_types, font_family, n_genes, gene_list=None):
    """
    Annotate the Plotly figure with cell type labels and horizontal bars.

    This function adds text annotations to the figure, representing the cell types. 
    If a gene list is not provided, it also adds horizontal bars with cell type 
    labels above the dot plot, providing a visual summary of the number of top genes 
    for each cell type.

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the annotations will be added.

    cell_types : list
        List of cell types to annotate.

    font_family : str
        Font family to use for the text annotations.

    n_genes : int
        Number of top genes considered. Used to determine the length of the horizontal bars.

    gene_list : list, optional
        List of specific genes considered. If provided, horizontal bars are not added.

    Returns:
    -------
    None
    """
    # Initialize variables for bar plotting
    bar_y_start = len(cell_types) + 0.1  # The y-position where the bars will start
    bar_x_start = 0  # The x-position where the first bar will start
    for idx, cell_type in enumerate(cell_types):
        fig.add_annotation(
            x=-1.25,
            y=idx,
            xref="x",
            yref="y",
            text=cell_type,
            showarrow=False,
            font=dict(
                size=12,
                color="black",
                family=font_family
            ),
            row=1,
            col=1
        )

    if gene_list is None:

        # Add the horizontal bars and vertical lines on top of the scatter plot
        for j, cell_type in enumerate(cell_types):
            bar_x_end = bar_x_start + n_genes - 1  # The x-position where the bar will end
            mid_x = (bar_x_start + bar_x_end) / 2  # Midpoint for text positioning

            # Add the horizontal bar
            fig.add_trace(
                go.Scatter(
                    x=[bar_x_start, bar_x_end],
                    y=[bar_y_start, bar_y_start],
                    mode='lines',
                    line=dict(color="black", width=1),
                    hoverinfo='none'
                ),
                row=1, col=1
            )

            # Add the cell type text at the midpoint of the bar
            fig.add_trace(
                go.Scatter(
                    x=[mid_x],
                    y=[bar_y_start],
                    mode='text',
                    text=[cell_type],
                    textposition="top center",
                    hoverinfo='none'
                ),
                row=1, col=1
            )

            # Add vertical lines (square brackets)
            for x_pos in [bar_x_start, bar_x_end]:
                fig.add_trace(
                    go.Scatter(
                        x=[x_pos, x_pos],
                        y=[bar_y_start - 0.2, bar_y_start + 0.02],
                        mode='lines',
                        line=dict(color="black", width=1),
                        hoverinfo='none'
                    ),
                    row=1, col=1
                )

            bar_x_start = bar_x_end + 1

    fig.add_annotation(
        x=-2.5,
        y=len(cell_types) / 2,
        text="Cell Types",
        textangle=-90,
        showarrow=False,
        font=dict(
            size=14,
            color="black",
            family=font_family
        ),
        row=1,
        col=1
    )
    fig.update_xaxes(fixedrange=True, domain=[0, 0.99], row=1, col=1)
    fig.update_yaxes(fixedrange=True, range=[-1, bar_y_start + 1], row=1, col=1)
    fig.update_xaxes(fixedrange=True, domain=[0.95, 0.985], row=1, col=2)
    fig.update_yaxes(fixedrange=True, range=[-5, 105], row=1, col=2)

def set_plot_layout(fig, all_genes, cell_types, font_family, base_height, base_width, scaling_factor):
    """
    Configure the layout and style of a Plotly dot plot figure.

    This function sets the appearance and layout of the figure, such as axis titles, 
    tick values, figure dimensions, and font styling. 

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure whose layout will be configured.

    all_genes : list
        List of all genes considered, used to set the x-axis tick values.

    cell_types : list
        List of cell types considered, used for layout scaling calculations.

    font_family : str
        Font family to be used for all text elements in the figure.

    base_height : int
        Base height of the figure.

    base_width : int
        Base width of the figure before scaling.

    scaling_factor : float
        Factor by which the figure width should be scaled, based on the number of genes.

    Returns:
    -------
    dict
        A dictionary of configuration options for interactive features of the figure.
    """
    adjusted_width = int(base_width * scaling_factor)
    config = {
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': 'custom_image',
            'height': base_height,
            'width': adjusted_width,
            'scale': 1
        }
    }
    fig.update_layout(
        xaxis1=dict(tickvals=list(range(len(all_genes))), ticktext=all_genes, title="Genes", showgrid=False, zeroline=False, tickangle=270),
        yaxis=dict(showgrid=False,zeroline=False,showline=False,ticks='',showticklabels=False),
        yaxis2=dict(showticklabels=False, showgrid=False, zeroline=False),
        xaxis2=dict(showticklabels=False, showgrid=False, zeroline=False),
        showlegend=False,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=1, r=30, t=20, b=50),
        font=dict(family=font_family, color="black"),
        width=adjusted_width,
        height=base_height
    )
    return config

def add_scatter_plot(fig, mean_expression_data, fraction_expressing_data, all_genes, cell_types, colorscale, max_size, flavor):
    """
    Add a scatter plot representation of gene expression data to a Plotly figure.

    This function generates and adds to the figure a scatter plot showing gene 
    expression across cell types. The size of each point represents the fraction 
    of cells expressing the gene, while the color indicates the mean expression 
    level or fold change, depending on the flavor.

    Parameters:
    ----------
    fig : plotly.graph_objects.Figure
        The Plotly figure to which the scatter plot will be added.

    mean_expression_data : numpy.ndarray
        2D array containing mean expression data of genes across cell types.

    fraction_expressing_data : numpy.ndarray
        2D array containing the fraction of cells expressing each gene across cell types.

    all_genes : list
        List of all genes considered, used for x-axis positioning.

    cell_types : list
        List of cell types, used for y-axis positioning.

    colorscale : str
        Colorscale for the scatter plot points, representing expression level or fold change.

    max_size : int
        Maximum size for the scatter plot points.

    flavor : str
        Specifies the method for obtaining gene values: 'mean_expression' or 'fold_change'.
        Determines the data visualized by the scatter plot's color.

    Returns:
    -------
    None
    """
    normalized_sizes = (fraction_expressing_data / np.max(fraction_expressing_data)) * max_size
    x_data, y_data, sizes, colors, hover_texts = [], [], [], [], []
    
    for i, gene in enumerate(all_genes):
        for j, cell_type in enumerate(cell_types):
            x_data.append(i)
            y_data.append(j)
            sizes.append(normalized_sizes[i, j])
            colors.append(mean_expression_data[i, j])
            hover_texts.append(f"Gene: {gene}<br>Fraction: {fraction_expressing_data[i, j]*100:.2f}%<br>Value: {mean_expression_data[i, j]:.2f}")
    
    final_colorscale = colorscale if flavor == "mean_expression" else "balance"
    color_range = [0, np.max(mean_expression_data)] if flavor == "mean_expression" else [-4, 4]
    
    fig.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            mode='markers',
            marker=dict(
                size=sizes,
                color=colors,
                colorscale=final_colorscale,
                showscale=False,
                colorbar=dict(title="", orientation='h'),
                cmin=color_range[0],
                cmax=color_range[1],
                line=dict(width=1, color="Black"),
                sizemode="diameter",
                opacity=1
            ),
            hoverinfo="text",
            hovertext=hover_texts
        ),
        row=1, col=1
    )


def plotly_dotplot_unified(adata, n_genes=10, gene_list=None, group_by="leiden", 
                           max_size=15, show_dendrogram=True, dendrogram_color='black', 
                           font_family='Arial', colorscale="Reds", flavor="mean_expression", min_logfc=1):
    """
    Create a unified dot plot visualization of gene expression using Plotly.

    This function produces a dot plot showcasing gene expression across different cell types.
    Each dot's size represents the fraction of cells expressing the gene, while its color 
    indicates the mean expression level or fold change. Additional annotations and an optional 
    dendrogram can also be added to provide context about the hierarchical structure of the data.

    Parameters:
    ----------
    adata : AnnData
        An object containing gene expression data, differential expression results, and 
        potentially dendrogram data.

    n_genes : int, optional
        Number of top genes to consider if gene_list is not provided. Default is 10.

    gene_list : list, optional
        List of specific genes to consider. If provided, n_genes is ignored.

    group_by : str, optional
        Column name in adata.obs that contains cell type or cluster information. Default is 'leiden'.

    max_size : int, optional
        Maximum size for the scatter plot points. Default is 15.

    show_dendrogram : bool, optional
        Whether to overlay a dendrogram on the plot, indicating hierarchical clustering. Default is True.

    dendrogram_color : str, optional
        Color for the dendrogram lines. Default is 'black'.

    font_family : str, optional
        Font family to use for all text elements in the figure. Default is 'Arial'.

    colorscale : str, optional
        Colorscale for the scatter plot points, representing expression level or fold change. Default is "Reds".

    flavor : str, optional
        Specifies the method for obtaining gene values: 'mean_expression' or 'fold_change'. Default is 'mean_expression'.

    min_logfc : float, optional
        Minimum log fold change to consider when selecting genes based on fold change. Only used if flavor is 'fold_change'. Default is 1.

    Returns:
    -------
    plotly.graph_objects.Figure
        The fully constructed Plotly figure.
    """
    
    cell_types = list(adata.uns["dendrogram_leiden"]['categories_ordered'])
    
    # Handling gene_list logic
    if gene_list is not None:
        n_genes = len(gene_list)  # Update n_genes to match the length of gene_list if provided
    
    mean_expression_data, fraction_expressing_data, all_genes = prepare_data(
        adata, cell_types, n_genes, gene_list, group_by, flavor, min_logfc)   # Pass min_logfc here
    
    fig = make_subplots(rows=1, cols=2, shared_yaxes=False, column_widths=[0.95, 0.05])
    
    # Add scatter plot
    add_scatter_plot(fig, mean_expression_data, fraction_expressing_data, all_genes, cell_types, colorscale, max_size, flavor)
    
    # Add Layout
    base_height = 350
    
    if gene_list is not None:
        if len(gene_list) < 20:
            base_width = 75
            scaling_factor = 5
        else:
            base_width = 55
            scaling_factor = 15
    else:
        scaling_factor = max(1, n_genes / 2)
        base_width = 500
    config = set_plot_layout(fig, all_genes, cell_types, font_family, base_height, base_width, scaling_factor)
    
    # Add Annotations
    add_annotations(fig, cell_types, font_family, n_genes, gene_list=gene_list)
    
    # Add dendrogram if required
    if show_dendrogram:
        add_dendrogram_to_plot(fig, adata, color_mapping={dendrogram_color: 'black'})
    
    fig.show(config=config)