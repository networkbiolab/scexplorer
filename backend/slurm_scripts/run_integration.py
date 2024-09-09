import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)

try:
    from utils import perform_batch_correction
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
def main(file_path, save_path,batch_correction_method, unique_id, preprocess):
    integration_plots_paths = perform_batch_correction(file_path, save_path,str(batch_correction_method),str(unique_id),bool(preprocess))
    print(f"{integration_plots_paths}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.argv[5])