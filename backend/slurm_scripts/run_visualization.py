import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)
try:
    from utils import genes_visualization
except ImportError as e:
    print(f"Importación fallida: {e}")
    
def main(file_path, save_path, dim_red, gene_list):
    list_gene = eval(gene_list)
    genes_plots_path = genes_visualization(file_path, save_path, dim_red, list_gene)
    print(f"{genes_plots_path}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])