import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)

try:
    from utils import load_and_preprocess_data
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
    
def main(actual_file_path, preprocess_plots_path,doublet_detection,root_path,unique_id, min_genes, min_cells,mito_threshold):
    preprocess_plots_paths = load_and_preprocess_data(actual_file_path, preprocess_plots_path,bool(doublet_detection),str(root_path),str(unique_id), int(min_genes), int(min_cells),int(mito_threshold))
    print(f"{preprocess_plots_paths}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])