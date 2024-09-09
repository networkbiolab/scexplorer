
import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)
try:
    from utils import get_quality_control_plots
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
def main(file_path, save_path, species, flavor,root_path,unique_id):
    plot_paths, dataset_description = get_quality_control_plots(file_path, save_path, species, flavor,root_path,str(unique_id))
    print(f"Plot Paths: {plot_paths}")
    print(f"Dataset Description: {dataset_description}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6])