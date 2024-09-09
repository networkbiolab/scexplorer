
import sys
from utils import get_quality_control_plots



def main(file_path, save_path, species, flavor):
    plot_paths, dataset_description = get_quality_control_plots(file_path, save_path, species, flavor)
    print(f"Plot Paths: {plot_paths}")
    print(f"Dataset Description: {dataset_description}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])