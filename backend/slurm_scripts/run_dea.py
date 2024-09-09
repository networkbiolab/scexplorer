import sys
import os
import json 

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)

try:
    from utils import dea
except ImportError as e:
    print(f"Importación fallida: {e}")
    
def main(actual_file_path, dea_plots_directory, n_genes, flavor, gene_list,method,root_path,unique_id):
    if gene_list:  # This will check if gene_list is not an empty string and not None
        list_gene = json.loads(gene_list)
    else:
        list_gene = ""

    dea_plots_paths = dea(actual_file_path, dea_plots_directory, int(n_genes), flavor, list_gene,str(method),root_path,unique_id)
    print(f"{dea_plots_paths}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])