import sys
import os

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
utils_path = os.path.join(parent_dir, 'utils.py')

sys.path.append(parent_dir)

try:
    from utils import var_names_adata
except ImportError as e:
    print(f"Importación fallida: {e}")
    
    
def main(file_path,):
    gene_list = var_names_adata(file_path)
    print(f"{gene_list}")

if __name__ == "__main__":
    main(sys.argv[1])