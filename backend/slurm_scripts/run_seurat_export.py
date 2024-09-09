import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)
try:
    from utils import run_r_script
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
def main(uuid):
    seurat_file_path =  run_r_script(uuid)
    print(f"{seurat_file_path}")


if __name__ == "__main__":
    main(sys.argv[1])