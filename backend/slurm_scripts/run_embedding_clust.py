import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)

try:
    from utils import run_clustering_on_data
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
def main(file_path, save_path, n_neighbors,n_pcs,resolution,root_path,unique_id):
    clustering_plots_paths = run_clustering_on_data(file_path, save_path, int(n_neighbors),n_pcs,float(resolution),root_path,unique_id)
    print(f"{clustering_plots_paths}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],sys.argv[6], sys.argv[7])