import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)

try:
    from utils import run_pca_on_data
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
def main(file_path, save_path, top_genes,hvg_flavor,root_path,unique_id):
    pca_plots_paths = run_pca_on_data(file_path, save_path, int(top_genes),hvg_flavor,root_path,unique_id)
    print(f"{pca_plots_paths}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])