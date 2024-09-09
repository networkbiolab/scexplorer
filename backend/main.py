import asyncio
import ast
import glob
import gzip
import json
import pandas as pd
import rpy2.robjects as robjects
import re
import os
import scanpy as sc
import shutil
import sqlite3
import subprocess
import uuid
import zipfile

from fastapi.staticfiles import StaticFiles
from fastapi import FastAPI, UploadFile, File, Depends, HTTPException,Form
from fastapi.responses import FileResponse,JSONResponse
from starlette.middleware.cors import CORSMiddleware
from pydantic import ValidationError
from utils import  upload_notification
from validations import ProcessDataParams,EmbeddingParams,PcaParams,VisualizationInput
from typing import Optional, List, Dict,Union




app = FastAPI(
    title="Single-Cell RNA-Seq Analysis API",
    description="""
    This API provides a comprehensive suite of tools for single-cell RNA sequencing data analysis.
    It offers functionalities including data upload, quality control, dimensionality reduction,
    clustering, differential expression analysis, and various visualization options.

    Key features:
    - File upload and processing
    - Quality control and filtering
    - PCA and embedding generation (UMAP)
    - Clustering analysis
    - Differential expression analysis
    - Gene expression visualization
    - Export options to Scanpy and Seurat formats
    - Comprehensive plotting and data retrieval endpoints

    For detailed usage instructions and parameter descriptions, please refer to the endpoint
    documentation below.
    """,
    version="1.0.0",
    terms_of_service="",
    contact={
        "name": "Sergio Hernández",
        "url": "https:tbd.cl",
        "email": "felipehg1991@gmail.com",
    },
    license_info={
        "name": "Apache 2.0",
        "url": "https://www.apache.org/licenses/LICENSE-2.0.html",
    },
    root_path="/backend"
)

UPLOAD_FOLDER = "uploads"
# os.makedirs(UPLOAD_FOLDER)
DATABASE = "tasks.db"


def get_db():
    conn = sqlite3.connect(DATABASE)
    return conn

@app.on_event("startup")
def startup():
    if not os.path.exists(UPLOAD_FOLDER):
        os.makedirs(UPLOAD_FOLDER)

    conn = get_db()
    cursor = conn.cursor()

    # Create tasks table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS tasks (
        id INTEGER PRIMARY KEY,
        uuid TEXT NOT NULL UNIQUE,
        file_path TEXT,
        species TEXT,
        email TEXT,
        analysisName TEXT,  
        created_at TEXT NOT NULL
    );
    """)

    # Create processing_info table with the new columns
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS processing_info (
        id INTEGER PRIMARY KEY,
        task_id INTEGER,
        cells_number_pre_processing INTEGER,
        genes_number_pre_processing INTEGER,
        cells_number_post_processing INTEGER,
        genes_number_post_processing INTEGER,
        min_genes INTEGER,
        min_cells INTEGER,
        mito_threshold REAL, 
        n_neighbors INTEGER DEFAULT NULL,
        n_pcs INTEGER DEFAULT NULL,
        clusters INTEGER DEFAULT NULL,
        FOREIGN KEY(task_id) REFERENCES tasks(id)
    );
    """)
        
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS dea_results (
        id INTEGER PRIMARY KEY,
        task_id INTEGER,
        dotplot_1_path TEXT,
        dotplot_2_path TEXT,
        FOREIGN KEY(task_id) REFERENCES tasks(id)
    );
    """)
    conn.commit()
    conn.close()
    
# app = FastAPI()
app.mount("/uploads", StaticFiles(directory="uploads"), name="uploads")


def normalize_files(dataset_3: Union[List[UploadFile], UploadFile, None] = File(None)) -> List[UploadFile]:
    if isinstance(dataset_3, list):
        return dataset_3
    elif isinstance(dataset_3, UploadFile):
        return [dataset_3]
    return [dataset_3]

def normalize_files(files: Union[List[UploadFile], UploadFile, None] = File(None)) -> List[UploadFile]:
    if isinstance(files, list):
        return files
    elif isinstance(files, UploadFile):
        return [files]
    return [files]

@app.get("/backend/")
def read_root():
    return {"Backend": "Running"}

@app.post("/backend/data_integration/", tags=["Main App"])
async def data_integration(
    dataset_1: List[UploadFile] = File(...),
    dataset_2: List[UploadFile] = File(...),
    dataset_3: List[UploadFile] = Depends(normalize_files),
    analysisName: str = Form(...),
    email: Optional[str] = Form(None),
    preprocess: bool = True,
    batch_correction_method: str = Form(...)
):
    """
    Upload files for single-cell RNA sequencing data integration.

    Parameters:
    - dataset_1: List[UploadFile] - First dataset files
    - dataset_2: List[UploadFile] - Second dataset files
    - dataset_3: Optional[List[UploadFile]] - Third dataset files (optional)
    - analysisName: str - Name of the analysis or dataset
    - batch_correction_method: str - Batch correction method to use ('combat', 'scanorama', 'bbknn', 'harmony')

    Returns:
    - A dictionary containing:
    - file_paths: List[str] - Paths where the files are stored
    - uuid: str - Unique identifier for the analysis
    - analysisName: str - Name of the analysis
    - plots_path: List[str] - Paths to generated integration plots (placeholder)
    - dataset_description: dict - Basic statistics about the dataset
    """
    try:
        unique_id = str(uuid.uuid4())
        dir_path = os.path.join(UPLOAD_FOLDER, unique_id)
        integration_plots_path = os.path.join(dir_path, "integration")

        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        if not os.path.exists(integration_plots_path):
            os.makedirs(integration_plots_path)

        file_paths = []
        datasets = [dataset_1, dataset_2] 
        print("dataset_3:", dataset_3)
        if dataset_3 != None:
            datasets.extend([dataset_3])

            
        print("datasets:", datasets)

        adata_list = []
        keys = []

        for i, dataset in enumerate(datasets):
            keys.append(f'batch{i+1}')
            for file in dataset:
                if file:
                    save_path = os.path.join(dir_path, file.filename)
                    with open(save_path, "wb+") as f:
                        f.write(file.file.read())
                    file_paths.append(save_path)

                    # Process files to AnnData objects
                    if save_path.endswith('.h5ad'):
                        adata = sc.read_h5ad(save_path)
                    elif save_path.endswith('.rds'):
                        r_script = f"""
                        library(Seurat)
                        library(zellkonverter)

                        # Specify the paths
                        seurat_object_path <- "{save_path}"
                        output_path <- "{save_path.replace('.rds', '.h5ad')}"

                        # Load the Seurat object
                        seurat_object <- readRDS(seurat_object_path)
                        seurat_object <- UpdateSeuratObject(seurat_object)

                        sce_object <- as.SingleCellExperiment(seurat_object)
                        writeH5AD(sce_object, output_path)

                        message("Seurat object has been successfully converted to H5AD format.")
                        """
                        r_script_path = os.path.join(dir_path, "convert_to_h5ad.R")
                        with open(r_script_path, "w") as f:
                            f.write(r_script)

                        cmd = ["Rscript", r_script_path]
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        if result.returncode != 0:
                            raise RuntimeError(result.stderr)

                        adata = sc.read_h5ad(save_path.replace('.rds', '.h5ad'))
                    else:
                        adata = sc.read_10x_mtx(
                            save_path,
                            var_names='gene_symbols',
                            cache=False
                        )
                    adata_list.append(adata)

        adata_concat = sc.concat(adata_list, label='batch', keys=keys)
        adata_concat.write_h5ad(os.path.join(dir_path, f"{unique_id}.h5ad"))
        file_path = os.path.join(dir_path, f"{unique_id}.h5ad")
    
        cmd = [
            "sbatch",
            "--job-name=integration_plots",
            f"--output={os.path.join(UPLOAD_FOLDER, unique_id, 'integration/integration_plots.out')}",
            f"--error={os.path.join(UPLOAD_FOLDER, unique_id, 'integration/integration_plots.err')}",
            "slurm_scripts/slurm_integration.sh",
            file_path,
            integration_plots_path,
            batch_correction_method,
            unique_id,
            str(preprocess)
        ]
        slurm_job_id = subprocess.run(cmd, capture_output=True, text=True).stdout.strip().split()[-1]

        job_status = await wait_for_job_completion(slurm_job_id)
        if job_status == "COMPLETED":
            if email:
                upload_notification(email, analysisName, unique_id)

            output_file_path = f"uploads/{unique_id}/integration/integration_plots.out"
            if os.path.exists(output_file_path):
                with open(output_file_path, 'r') as file:
                    output_contents = file.read()

        else:
            raise HTTPException(status_code=404, detail="Output file not found.")   
        
        file_paths_list = re.findall(r"\['[^]]*'\]", output_contents)
        if file_paths_list:
            file_paths_list = eval(file_paths_list[-1])
            print(f"file_paths_list: {file_paths_list}")
        else:
            file_paths_list = []
        return {
            "message": "File uploaded successfully",
            "uuid": unique_id,
            "analysisName": analysisName,
            "plots_path": file_paths_list,
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    

@app.exception_handler(Exception)
async def exception_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={"message": "An unexpected error occurred", "details": str(exc)},
    )

def run_r_script(script_path, script_args):
    cmd = ["Rscript", script_path] + script_args
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return result.stdout

@app.post("/backend/upload/", tags=["Main App"])
async def upload_file(
    species: str = Form(...),
    analysisName: str = Form(...),
    email: Optional[str] = Form(None),
    flavor: str = Form(...),
    files: List[UploadFile] = File(None)
):
    """
    Upload files for single-cell RNA sequencing analysis.

    Parameters:
    - files: List[UploadFile] - The single-cell RNA-seq data files to be uploaded
    - species: str - The species of the sample (e.g., "human", "mouse")
    - email: Optional[str] - Email address for notifications (optional)
    - analysisName: str - Name of the analysis or dataset

    Returns:
    - A dictionary containing:
    - file_paths: List[str] - Paths where the files are stored
    - uuid: str - Unique identifier for the analysis
    - analysisName: str - Name of the analysis
    - plots_path: List[str] - Paths to generated QC plots
    - dataset_description: dict - Basic statistics about the dataset
    """
    try:
        example_files = {
            "example_human": "examples/pbmc35.h5ad",
            "example_mouse": "examples/mouse_brain.h5ad",
            "example_zebrafish": "examples/zebrafish.h5ad",
        }

        if analysisName.startswith("example_"):
            unique_id = str(uuid.uuid4())
            dir_path = os.path.join(UPLOAD_FOLDER, unique_id)
            example_path = analysisName
            file_path = example_files.get(example_path)
            if not file_path:
                raise HTTPException(status_code=400, detail="Example not found")
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
            qc_plots_path = os.path.join(dir_path, "qc_plots")
            if not os.path.exists(qc_plots_path):
                os.makedirs(qc_plots_path)
            new_filename = f"{unique_id}.h5ad"
            new_file_path = os.path.join(dir_path, new_filename)
            shutil.copy(file_path, new_file_path)
        else:
            # If not an example, handle new file uploads
            unique_id = str(uuid.uuid4())
            dir_path = os.path.join(UPLOAD_FOLDER, unique_id)
            
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
            qc_plots_path = os.path.join(dir_path, "qc_plots")
            if not os.path.exists(qc_plots_path):
                os.makedirs(qc_plots_path)
            
            file_paths = []
            features_path = None

            if len(files) == 1:
                file = files[0]
                file_extension = os.path.splitext(file.filename)[-1]
                save_path = os.path.join(dir_path, f"{unique_id}{file_extension}")
                with open(save_path, "wb+") as f:
                    f.write(file.file.read())
                file_paths.append(save_path)
            else:
                for file in files:
                    print(f"filename:{file.filename}")
                    save_path = os.path.join(dir_path, file.filename)
                    with open(save_path, "wb+") as f:
                        f.write(file.file.read())
                    file_paths.append(save_path)
                    if 'features.tsv.gz' in file.filename:
                        print(f"features detected:{file.filename}")
                        features_path = save_path

            if any(file.endswith('.h5ad') for file in file_paths):
                file_path = next(file for file in file_paths if file.endswith('.h5ad'))
            elif any(file.endswith('.h5') for file in file_paths):
                print("inside h5")
                print(f"file: {file_paths}")
                try:
                    adata = sc.read_10x_h5(file_paths[0])
                    file_path = os.path.join(dir_path, f"{unique_id}.h5ad")
                    adata.write_h5ad(file_path)
                except Exception as e:
                    print(f"error: {e}")
            elif any(file.endswith('.rds') for file in file_paths):
                file_path = next(file for file in file_paths if file.endswith('.rds'))

                r_script = f"""
                library(Seurat)
                library(zellkonverter)

                # Specify the paths
                seurat_object_path <- "{os.path.join(dir_path, unique_id)}.rds"
                output_path <- "{os.path.join(dir_path, unique_id)}.h5ad"

                # Load the Seurat object
                seurat_object <- readRDS(seurat_object_path)
                seurat_object <- UpdateSeuratObject(seurat_object)

                sce_object <- as.SingleCellExperiment(seurat_object)
                writeH5AD(sce_object, output_path)

                message("Seurat object has been successfully converted to H5AD format.")
                """
                r_script_path = os.path.join(dir_path, "convert_to_h5ad.R")
                with open(r_script_path, "w") as f:
                    f.write(r_script)
                run_r_script(r_script_path, [])
                file_path = os.path.join(dir_path, f"{unique_id}.h5ad")
            else:
                print("inside_10x")
                print(f"dir_path:{dir_path}")

                try:
                    if features_path:
                        print(f"inside {features_path}")
                        with gzip.open(features_path, 'rt') as f:
                            features_df = pd.read_csv(f, sep='\t', header=None)
            
                        # Check if features_df has exactly one column
                        if features_df.shape[1] == 1:
                            features_df[1] = features_df[0]  # Copy the first column into the second
                            features_df[2] = "Gene Expression"  # Create a third column named "Gene Expression"
                        elif features_df.shape[1] == 2:
                            features_df[2] = "Gene Expression"  # Add a third column named "Gene Expression"
                        
                        # Write the updated DataFrame back to a new file and recompress it
                        new_features_path = features_path.replace('.gz', '')
                        features_df.to_csv(new_features_path, sep='\t', header=False, index=False)
                        with open(new_features_path, 'rb') as f_in:
                            with gzip.open(features_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(new_features_path)

                        
                        if os.path.exists(qc_plots_path) and os.path.isdir(qc_plots_path):
                            shutil.rmtree(qc_plots_path)
                    
                    print(f"files:{os.listdir(dir_path)}")
                    
                    adata = sc.read_10x_mtx(
                    dir_path,  # directory with matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
                    var_names='gene_symbols',  # use gene symbols for variable names (otherwise use gene IDs)
                    )
                    if not os.path.exists(qc_plots_path):
                        os.makedirs(qc_plots_path)
                    print(f"adata: {adata}")
                    adata.write_h5ad(os.path.join(dir_path, f"{unique_id}.h5ad"))
                    print("post_h5ad_10x writing")
                    file_path = os.path.join(dir_path, f"{unique_id}.h5ad")
                except Exception as e:
                    print(f"error:{e}")
        print("pre_slurm")
        root_path = UPLOAD_FOLDER
        cmd = [
            "sbatch",
            "--job-name=qc_plots",
            f"--output={os.path.join(UPLOAD_FOLDER, unique_id, 'qc_plots/qc_plots.out')}",
            f"--error={os.path.join(UPLOAD_FOLDER, unique_id, 'qc_plots/qc_plots.err')}",
            "slurm_scripts/qc_plots.sh",
            file_path,
            qc_plots_path,
            species,
            flavor,
            root_path,
            unique_id
        ]
        slurm_job_id = subprocess.run(cmd, capture_output=True, text=True).stdout.strip().split()[-1]

        job_status = await wait_for_job_completion(slurm_job_id)
        if job_status == "COMPLETED":
            if email:
                upload_notification(email, analysisName, unique_id)

            output_file_path = os.path.join(UPLOAD_FOLDER, unique_id, "qc_plots/qc_plots.out")
            print("post_slurm")

            if os.path.exists(output_file_path):
                with open(output_file_path, 'r') as f:
                    output_contents = f.read()
                plot_paths, dataset_description = parse_output_content(output_contents)
            else:
                raise HTTPException(status_code=404, detail="Output file not found.")
            # if file_paths ==None:
            #     file_paths = []
            return {
                "message": "File uploaded successfully",
                # "file_paths": file_paths,
                "uuid": unique_id,
                "analysisName": analysisName,
                "plots_path": plot_paths,
                "dataset_description": dataset_description
            }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
@app.get("/backend/process/{file_uuid}/", tags=["Main App"])
async def process_file(file_uuid: str, params: ProcessDataParams = Depends()):
    """
    Process a previously uploaded single-cell RNA-seq file with quality control filters.

    Parameters:
    - file_uuid: str - UUID of the uploaded file
    - params: ProcessDataParams
    - min_genes: int - Minimum number of genes expressed for a cell to be included (default: 200)
    - min_cells: int - Minimum number of cells a gene must be expressed in to be included (default: 3)
    - mito_threshold: int - Maximum percentage of mitochondrial genes allowed (default: 5)

    Returns:
    - A dictionary containing processing results, including:
    - Cells and genes counts before and after filtering
    - Paths to QC plots
    """
    dir_path = os.path.join(UPLOAD_FOLDER, file_uuid)
    file_path = os.path.join(dir_path, f"{file_uuid}.h5ad")
    matching_files = glob.glob(file_path)
    
    if not matching_files:
        return {"error": "Archivo no encontrado."}
    
    preprocess_plots_path = os.path.join(dir_path, "preprocess")
        
    if not os.path.exists(preprocess_plots_path):
        os.makedirs(preprocess_plots_path)
        
    actual_file_path = matching_files[0]
    
    partition = "debug" 
    cpus_per_task = 16  
    memory = "16G"  
    root_path = UPLOAD_FOLDER
    cmd = [
        "sbatch",
        "--job-name=preprocess",
        f"--output={os.path.join(UPLOAD_FOLDER, file_uuid, 'preprocess/preprocess.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, file_uuid, 'preprocess/preprocess.err')}",
        "--ntasks=1",
        f"--cpus-per-task={cpus_per_task}",  
        f"--mem={memory}",  
        f"--partition={partition}",  
        "slurm_scripts/slurm_preprocess.sh",
        actual_file_path,
        preprocess_plots_path,
        str(params.doublet_detection),
        str(root_path),
        str(file_uuid),
        str(params.min_genes),
        str(params.min_cells),
        str(params.mito_threshold)
 
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise HTTPException(status_code=500, detail=f"Error submitting SLURM job: {result.stderr}")
    
    # Debugging output to check the result of the sbatch command
    print(f"sbatch result: {result.stdout.strip()}")
    
    # Extracting job ID correctly
    slurm_job_id_line = result.stdout.strip().split()[-1]
    try:
        slurm_job_id = int(slurm_job_id_line)
    except ValueError:
        raise HTTPException(status_code=500, detail=f"Failed to extract job ID from sbatch output: {slurm_job_id_line}")
    
    print(f"SLURM JOB ID: {slurm_job_id}")  # Print the job ID for debugging
    
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""

    if job_status == "COMPLETED":
        output_file_path = f"uploads/{file_uuid}/preprocess/preprocess.out"
        output_contents = await read_output_file(output_file_path)
    else:
        raise HTTPException(status_code=500, detail=f"Job did not complete successfully. Status: {job_status}")
    
    match = re.search(r'({.*})', output_contents, re.DOTALL)
    if match:
        dict_part = match.group(1)
        try:
            processing_results = eval(dict_part)
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Failed to evaluate dictionary: {e}")
    else:
        raise HTTPException(status_code=500, detail="Failed to extract processing results from output.")
    
    return processing_results


async def read_output_file(output_file_path):
    retries = 0
    max_retries = 10
    while retries < max_retries:
        if os.path.exists(output_file_path):
            with open(output_file_path, 'r') as file:
                return file.read()
        retries += 1
        await asyncio.sleep(5)  # Espera 5 segundos antes de volver a verificar
    raise HTTPException(status_code=500, detail="Output file not found after job completion.")

@app.post("/backend/pca",tags=["Main App"])
async def pca(file_uuid: str,params : PcaParams = Depends()):
    """
    Perform Principal Component Analysis (PCA) on a processed single-cell RNA-seq dataset.

    Parameters:
    - file_uuid: str - UUID of the processed file
    - params: PcaParams
    - n_genes: int - Number of highly variable genes to use for PCA (default: 2000)
    - flavor: str - Method for identifying highly variable genes (e.g., "seurat", "cell_ranger")

    Returns:
    - A dictionary containing:
    - file_uuid: str - UUID of the processed file
    - pca_plot_paths: List[str] - Paths to PCA plots (e.g., elbow plot, PCA components plot)
    """
    dir_path = os.path.join(UPLOAD_FOLDER, file_uuid)
    file_path = os.path.join(dir_path, f"{file_uuid}.h5ad")
    matching_files = glob.glob(file_path)
    
    if not matching_files:
        return {"error": "Archivo no encontrado."}
    
    embedding_plots_directory = os.path.join(dir_path,"embedding")
    if not os.path.exists(embedding_plots_directory):
        os.makedirs(embedding_plots_directory)
        
    actual_file_path = matching_files[0]
    
    if not os.path.exists(actual_file_path):
        raise HTTPException(status_code=404, detail="File not found in the upload folder")
    
    partition = "debug" 
    cpus_per_task = 8  
    memory = "16G"
    root_path = UPLOAD_FOLDER
    print(f"backend pca: {actual_file_path}")
    cmd = [
        "sbatch",
        "--job-name=embedding",
        f"--output={os.path.join(UPLOAD_FOLDER, file_uuid, f'embedding/embedding.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, file_uuid, f'embedding/embedding.err')}",
        "--ntasks=1",
        f"--cpus-per-task={cpus_per_task}",  
        f"--mem={memory}",  
        f"--partition={partition}",
        "slurm_scripts/slurm_embedding_pca.sh",
        actual_file_path,
        embedding_plots_directory,
        str(params.n_genes),
        params.flavor,
        root_path,
        file_uuid
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"sbatch result: {result.stdout.strip()}")
    
    # Extracting job ID correctly
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}")
        
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""

    if job_status == "COMPLETED":

        output_file_path = f"uploads/{file_uuid}/embedding/embedding.out"
        output_contents = await read_output_file(output_file_path)

    else:
        raise HTTPException(status_code=404, detail="Output file not found.")   
    
    try:
        file_paths_list = eval(output_contents)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to return pca file paths: {e}")

    return {
    "file_uuid": file_uuid, 
    "pca_plot_paths": file_paths_list,
    "message": "PCA Analysis done."
    }
    
@app.post("/backend/embedding/",tags=["Main App"])
async def embedding(file_uuid: str,params: EmbeddingParams = Depends()): 
    """
    Perform dimensionality reduction and clustering on a processed single-cell RNA-seq dataset.

    Parameters:
    - file_uuid: str - UUID of the processed file
    - params: EmbeddingParams
    - n_neighbors: int - Number of neighbors for UMAP (default: 15)
    - n_pcs: int - Number of principal components to use (default: 18)

    Returns:
    - A dictionary containing:
    - file_uuid: str - UUID of the processed file
    - embedding_plot_paths: List[str] - Paths to embedding plots (e.g., UMAP, t-SNE)
    """
    dir_path = os.path.join(UPLOAD_FOLDER, file_uuid)
    file_path = os.path.join(dir_path, f"{file_uuid}.h5ad")
    matching_files = glob.glob(file_path)
    
    if not matching_files:
        return {"error": "Archivo no encontrado."}
    
    embedding_plots_directory = os.path.join(dir_path,"embedding")
    if not os.path.exists(embedding_plots_directory):
        os.makedirs(embedding_plots_directory)
        
    actual_file_path = matching_files[0]
    
    if not os.path.exists(actual_file_path):
        raise HTTPException(status_code=404, detail="File not found in the upload folder")
    
    root_path = UPLOAD_FOLDER

    cmd = [
        "sbatch",
        "--job-name=embedding_clustering",
        f"--output={os.path.join(UPLOAD_FOLDER, file_uuid, f'embedding/embedding_clustering.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, file_uuid, f'embedding/embedding_clustering.err')}",
        "slurm_scripts/slurm_embedding_clust.sh",
        actual_file_path,
        embedding_plots_directory,
        str(params.n_neighbors),
        str(params.n_pcs),
        str(params.resolution),
        root_path,
        file_uuid
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"sbatch result: {result.stdout.strip()}")
    
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}")
        
    job_status = await wait_for_job_completion(slurm_job_id)   
    output_contents = ""
     
    if job_status == "COMPLETED":
        output_file_path = f"uploads/{file_uuid}/embedding/embedding_clustering.out"
        output_contents = await read_output_file(output_file_path)

    else:
        raise HTTPException(status_code=404, detail="Output file not found.")   
    
    try:
        embedding_plot_paths = eval(output_contents)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to return umap file paths: {e}")
    
    return {
        "file_uuid": file_uuid, 
        "embedding_plot_paths": embedding_plot_paths,
        "message": "Embedding and clustering done."
    }
    
    # clusters,embedding_plot_paths = run_clustering_on_data(file_path=actual_file_path,n_neighbors=params.n_neighbors,n_pcs=params.n_pcs,save_path=embedding_plots_directory)
    
    # conn = get_db()
    # cursor = conn.cursor()
    # # Update the database
    # cursor.execute("""
    # UPDATE processing_info 
    # SET n_neighbors = ?, n_pcs = ?, clusters = ?
    # WHERE task_id = (SELECT id FROM tasks WHERE uuid = ?)
    # """, (params.n_neighbors, params.n_pcs, clusters, file_uuid))
    # conn.commit()
    # conn.close()
    
    # return {
    #     "file_uuid": file_uuid, 
    #     "clusters": clusters, 
    #     "embedding_plot_paths": embedding_plot_paths,
    #     "message": "Embedding and clustering done."
    # }
    
@app.post("/backend/genes_visualization/",tags=["Main App"])
async def visualization(vis:VisualizationInput):
    """
    Visualize gene expression for specified genes in a dimensionality-reduced space.

    Parameters:
    - vis: VisualizationInput
    - file_uuid: str - UUID of the processed file
    - dim_red: str - Dimensionality reduction method to use (e.g., "umap")
    - gene_list: List[str] - List of gene names to visualize

    Returns:
    - A dictionary containing paths to gene expression visualization plots
    """
    dir_path = os.path.join(UPLOAD_FOLDER, vis.file_uuid)
    file_path = os.path.join(dir_path, f"{vis.file_uuid}.h5ad")
    matching_files = glob.glob(file_path)
    
    if not matching_files:
        return {"error": "Archivo no encontrado."}
    
    genes_plots_directory = os.path.join(dir_path,"visualization")
    
    if not os.path.exists(genes_plots_directory):
        os.makedirs(genes_plots_directory)
    actual_file_path = matching_files[0]
    
    
    cmd = [
        "sbatch",
        "--job-name=visualization",
        f"--output={os.path.join(UPLOAD_FOLDER, vis.file_uuid, f'visualization/visualization.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, vis.file_uuid, f'visualization/visualization.err')}",
        "slurm_scripts/slurm_visualization.sh",
        actual_file_path,
        genes_plots_directory,
        vis.dim_red,
        str(vis.gene_list)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"sbatch result: {result.stdout.strip()}")
    
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}") 
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""

    if job_status == "COMPLETED":
        output_file_path = f"uploads/{vis.file_uuid}/visualization/visualization.out"
        output_contents = await read_output_file(output_file_path)

    else:
        raise HTTPException(status_code=404, detail="Output file not found.")   
    
    try:
        genes_plot_paths = eval(output_contents)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to return visualization file paths: {e}")
    
    return genes_plot_paths    
    
@app.post("/backend/gene_names",tags=["Main App"])
async def gene_names(uuid: str): 
    """
    Retrieve the list of gene names present in a processed single-cell RNA-seq dataset.

    Parameters:
    - uuid: str - UUID of the processed file

    Returns:
    - A list of gene names (strings) present in the dataset
    """
    dir_path = os.path.join(UPLOAD_FOLDER, uuid)
    file_path = os.path.join(dir_path, f"{uuid}.h5ad")
    matching_files = glob.glob(file_path)
    
    if not matching_files:
        return {"error": "Archivo no encontrado."}
    
    embedding_plots_directory = os.path.join(dir_path,"embedding")
    if not os.path.exists(embedding_plots_directory):
        os.makedirs(embedding_plots_directory)
        
    actual_file_path = matching_files[0]
    
    if not os.path.exists(actual_file_path):
        raise HTTPException(status_code=404, detail="File not found in the upload folder")
    
    
    cmd = [
        "sbatch",
        "--job-name=embedding_gene_names",
        f"--output={os.path.join(UPLOAD_FOLDER, uuid, f'embedding/embedding_gene_names.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, uuid, f'embedding/embedding_gene_names.err')}",
        "slurm_scripts/slurm_gene.sh",
        actual_file_path,
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"sbatch result: {result.stdout.strip()}")
    
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}") 
            
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""

    if job_status == "COMPLETED":
        output_file_path = f"uploads/{uuid}/embedding/embedding_gene_names.out"
        output_contents = await read_output_file(output_file_path)

    else:
        raise HTTPException(status_code=404, detail="Output file not found.")   
    
    try:
        output_contents = eval(output_contents)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to return genes names: {e}")
    
    return output_contents

@app.get("/backend/info/{file_uuid}/",tags=["Main App"])
async def get_info(file_uuid: str):
    """
    Retrieve detailed information about a specific single-cell RNA-seq analysis task.

    Parameters:
    - file_uuid: str - UUID of the analysis task

    Returns:
    - A dictionary containing task and processing information, including:
    - File path
    - Species
    - Email (if provided)
    - Analysis name
    - Creation date
    - Processing parameters (e.g., min_genes, min_cells, mito_threshold)
    - Cell and gene counts before and after processing
    - Clustering information (if available)
    """
    conn = get_db()
    cursor = conn.cursor()
    cursor.execute("""
    SELECT t.*, p.*
    FROM tasks t
    LEFT JOIN processing_info p ON t.id = p.task_id
    WHERE t.uuid = ?
    """, (file_uuid,))
    data = cursor.fetchone()

    # Get column names
    columns = [description[0] for description in cursor.description]
    conn.close()

    if not data:
        return {"error": "No data found for given UUID."}

    # Convert the tuple data to a dictionary with column names as keys
    result = dict(zip(columns, data))
    return result


@app.post("/backend/DEA/", tags=["Main App"])
async def DEA(file_uuid: str, n_genes: int, flavor: str, method : str = 'wilcoxon',gene_list: list[str] = None):
    """
    Perform Differential Expression Analysis on a processed single-cell RNA-seq dataset.

    Parameters:
    - file_uuid: str - UUID of the processed file
    - n_genes: int - Number of top differentially expressed genes to consider
    - flavor: str - Method for differential expression analysis (e.g., "t-test", "wilcoxon")
    - gene_list: Optional[List[str]] - List of specific genes to analyze (if None, uses top n_genes)

    Returns:
    - A dictionary containing:
    - Paths to DEA visualization plots (e.g., volcano plots, heatmaps)
    - Differential expression statistics
    """
    # print(f"Received: {file_uuid}, {n_genes}, {flavor}, {gene_list}") 
    dir_path = os.path.join(UPLOAD_FOLDER, file_uuid)
    file_path = os.path.join(dir_path, f"{file_uuid}.h5ad")
    matching_files = glob.glob(file_path)
    
    if not matching_files:
        return {"error": "File not found."}
    
    dea_plots_directory = os.path.join(dir_path, "dea")
    
    if not os.path.exists(dea_plots_directory):
        os.makedirs(dea_plots_directory)
    actual_file_path = matching_files[0]
    
    gene_list_str = json.dumps(gene_list) if gene_list else "[]"
    root_path = UPLOAD_FOLDER

    cmd = [
        "sbatch",
        "--job-name=dea",
        f"--output={os.path.join(UPLOAD_FOLDER, file_uuid, 'dea/dea.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, file_uuid, 'dea/dea.err')}",
        "slurm_scripts/slurm_dea.sh",
        actual_file_path,
        dea_plots_directory,
        str(n_genes),
        flavor,
        gene_list_str,
        method,
        root_path,
        file_uuid
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"sbatch result: {result.stdout.strip()}")
    
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}") 
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""
    output_file_path = f"uploads/{file_uuid}/dea/dea.out"
    if job_status == "COMPLETED" and os.path.exists(output_file_path):
        output_contents = await read_output_file(output_file_path)
        try:
            dea_plot_paths = eval(output_contents)
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Failed to return DEA path: {e}")
        return dea_plot_paths
    else:
        detail_msg = "Job did not complete successfully or output file was not found."
        raise HTTPException(status_code=404, detail=detail_msg)


@app.post("/backend/scanpy_export",tags=["Main App"])
async def export_to_scanpy(uuid:str):
    """
    Export the processed single-cell RNA-seq analysis results to Scanpy (.h5ad) format.

    Parameters:
    - uuid: str - UUID of the analysis task

    Returns:
    - A FileResponse containing the exported Scanpy (.h5ad) file
    """
    scanpy_path = f"uploads/{uuid}/{uuid}.h5ad" 
    if os.path.exists(scanpy_path):
        return FileResponse(path=scanpy_path, filename=f"{uuid}.h5ad")
    else:
        raise HTTPException(status_code=404, detail="File not found.") 
    
@app.post("/backend/seurat_export",tags=["Main App"])
async def export_to_seurat(uuid:str):
    """
    Export the processed single-cell RNA-seq analysis results to Seurat (.rds) format.

    Parameters:
    - uuid: str - UUID of the analysis task

    Returns:
    - A FileResponse containing the exported Seurat (.rds) file
    """
    
    
    cmd = [
        "sbatch",
        "--job-name=seurat_export",
        f"--output={os.path.join(UPLOAD_FOLDER, uuid, f'seurat_export.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, uuid, f'seurat_export.err')}",
        "slurm_scripts/slurm_seurat_export.sh",
        uuid
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"sbatch result: {result.stdout.strip()}")
    
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}") 
            
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""
    if job_status == "COMPLETED":
        output_file_path = f"uploads/{uuid}/seurat_export.out"
        output_contents = await read_output_file(output_file_path)

    else:
        raise HTTPException(status_code=404, detail="Output file not found.")   

    match = re.search(r"(uploads/[^\s]+)", output_contents)

    # Extracting the file path if present
    file_path = match.group(0) if match else None

    # Print the extracted file path

    if os.path.exists(file_path):
        return FileResponse(path=file_path, filename=f"{uuid}.rds")
    else:
        raise HTTPException(status_code=404, detail="File not found.") 
    
@app.post("/backend/clustree",tags=["Main App"])
async def clustree_R(uuid:str):
    """
    Generate a Clustree plot to visualize clustering stability across different resolutions.

    Parameters:
    - uuid: str - UUID of the analysis task

    Returns:
    - A FileResponse containing the Clustree plot (SVG format)
    """

    dir_path = f"uploads/{uuid}/"
    dea_plots_directory = os.path.join(dir_path,"embedding")
    
    if not os.path.exists(dea_plots_directory):
        os.makedirs(dea_plots_directory)
        
    cmd = [
        "sbatch",
        "--job-name=clustree",
        f"--output={os.path.join(UPLOAD_FOLDER, uuid, f'embedding/clustree.out')}",
        f"--error={os.path.join(UPLOAD_FOLDER, uuid, f'embedding/clustree.err')}",
        "slurm_scripts/slurm_seurat_clustree.sh",
        uuid
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True) 
    print(f"sbatch result: {result.stdout.strip()}")
    
    slurm_job_id = result.stdout.strip().split()[-1]
    print(f"SLURM JOB ID: {slurm_job_id}") 
            
    job_status = await wait_for_job_completion(slurm_job_id)
    output_contents = ""
    if job_status == "COMPLETED":
        output_file_path = f"uploads/{uuid}/embedding/clustree.out"
        output_contents = await read_output_file(output_file_path)


    else:
        raise HTTPException(status_code=404, detail="Output file not found.")   

    match = re.search(r"(uploads/[^\s]+)", output_contents)

    # Extracting the file path if present
    file_path = match.group(0) if match else None

    if os.path.exists(file_path):
        return FileResponse(path=file_path, filename="clustreePlot.svg")
    else:
        raise HTTPException(status_code=404, detail="File not found.") 
    
    
    
@app.get("/backend/export-files/{uuid}",tags=["Main App"])
async def export_files(uuid: str):
    """
    Export all analysis files and results for a single-cell RNA-seq analysis task as a ZIP archive.

    Parameters:
    - uuid: str - UUID of the analysis task

    Returns:
    - A FileResponse containing a ZIP file with all analysis plots, data, and results
    """
    base_path = f'uploads/{uuid}'
    if not os.path.exists(base_path):
        raise HTTPException(status_code=404, detail="UUID not found")

    # Define the name of the ZIP file
    zip_filename = f"/tmp/{uuid}.zip"
    
    # Create a ZIP file
    with zipfile.ZipFile(zip_filename, 'w') as zipf:
        # Walk through the directory structure
        for root, dirs, files in os.walk(base_path):
            for file in files:
                if file.endswith('.html') or file.endswith('.svg'):
                    # Create a complete file path
                    file_path = os.path.join(root, file)
                    # Write file to the zip archive
                    zipf.write(file_path, arcname=os.path.relpath(file_path, base_path))

    return FileResponse(zip_filename, media_type='application/zip', filename=f"{uuid}.zip")    
########################################################################     
@app.get("/backend/uploads/{plot_path:path}", tags=["Main App"])
async def serve_plot_raw_data(plot_path: str):
    full_path = f"/app/uploads/{plot_path}"
    # print(f"Trying to access: {full_path}")
    # print(f"File exists: {os.path.exists(full_path)}")
    # print(f"Directory contents: {os.listdir(os.path.dirname(full_path))}")
    # full_path =  plot_path
    # print(f"Trying to access: {full_path}")  

    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)


@app.get("/backend/plots_preprocessing/{plot_path:path}",tags=["Main App"])
async def serve_plot_plots_preprocessing(plot_path: str):
    """
    Serve raw plot data for a specific analysis task.

    Parameters:
    - plot_path: str - Path to the requested plot file

    Returns:
    - A FileResponse containing the requested plot file
    """
    full_path =  plot_path
    print(f"Trying to access: {full_path}") 
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)


@app.get("/backend/plots_pca_embedding/{plot_path:path}",tags=["Main App"])
async def serve_plot_plots_embedding(plot_path: str):
    """
    Serve raw plot data for a specific analysis task.

    Parameters:
    - plot_path: str - Path to the requested plot file

    Returns:
    - A FileResponse containing the requested plot file
    """
    full_path =  plot_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)

@app.get("/backend/plots_embedding/{plot_path:path}",tags=["Main App"])
async def serve_plot_plots_embedding(plot_path: str):
    """
    Serve raw plot data for a specific analysis task.

    Parameters:
    - plot_path: str - Path to the requested plot file

    Returns:
    - A FileResponse containing the requested plot file
    """
    full_path =  plot_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)

@app.get("/backend/plots_dea/{plot_path:path}",tags=["Main App"])
async def serve_plot_plots_dea(plot_path: str):
    """
    Serve raw plot data for a specific analysis task.

    Parameters:
    - plot_path: str - Path to the requested plot file

    Returns:
    - A FileResponse containing the requested plot file
    """
    full_path =  plot_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)

@app.get("/backend/genes_vis/{plot_path:path}",tags=["Main App"])
async def serve_plot_plot_clustree(plot_path: str):
    """
    Serve raw plot data for a specific analysis task.

    Parameters:
    - plot_path: str - Path to the requested plot file

    Returns:
    - A FileResponse containing the requested plot file
    """
    full_path =  plot_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)

@app.get("/backend/plots_clustree/{plot_path:path}",tags=["Main App"])
async def serve_plot_plot_clustree(plot_path: str):
    """
    Serve raw plot data for a specific analysis task.

    Parameters:
    - plot_path: str - Path to the requested plot file

    Returns:
    - A FileResponse containing the requested plot file
    """
    full_path =  plot_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="Plot not found.")
    return FileResponse(full_path)


@app.get("/backend/get_obs/{csv_path:path}", tags=["Main App"])
async def serve_csv_file(csv_path: str):
    """
    Serve raw CSV file for a specific analysis task.

    Parameters:
    - csv_path: str - Path to the requested CSV file

    Returns:
    - A FileResponse containing the requested CSV file
    """
    full_path = csv_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="CSV not found.")
    return FileResponse(full_path)


@app.get("/backend/get_var/{csv_path:path}", tags=["Main App"])
async def serve_csv_file(csv_path: str):
    """
    Serve raw CSV file for a specific analysis task.

    Parameters:
    - csv_path: str - Path to the requested CSV file

    Returns:
    - A FileResponse containing the requested CSV file
    """
    full_path = csv_path
    print(f"Trying to access: {full_path}")  
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="CSV not found.")
    return FileResponse(full_path)
##mover a utils despues



@app.get("/backend/list_csv_files/{uuid}", tags=["Main App"])
async def list_csv_files(uuid: str):
    """
    List all CSV files for a specific analysis task.

    Parameters:
    - uuid: str - UUID for the specific analysis

    Returns:
    - A list of filenames present in the directory
    """
    directory = os.path.join(UPLOAD_FOLDER, uuid, "dea", "results")
    if not os.path.exists(directory):
        raise HTTPException(status_code=404, detail="Directory not found.")
    
    files = [f for f in os.listdir(directory) if f.endswith('.csv')]
    return {"files": files}

@app.get("/backend/get_csv/{uuid}/{filename}", tags=["Main App"])
async def serve_csv_file(uuid: str, filename: str):
    """
    Serve raw CSV file for a specific analysis task.

    Parameters:
    - uuid: str - UUID for the specific analysis
    - filename: str - Name of the CSV file

    Returns:
    - A FileResponse containing the requested CSV file
    """
    directory = os.path.join(UPLOAD_FOLDER, uuid, "dea", "results")
    full_path = os.path.join(directory, filename)
    if not os.path.exists(full_path):
        raise HTTPException(status_code=404, detail="CSV not found.")
    return FileResponse(full_path)



def parse_scontrol_job_state(job_id):
    """Retrieve and parse the job state from scontrol."""
    cmd = ["scontrol", "show", "jobid", f"{str(job_id)}"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    # print(f"scontrol : {result}")
    if result.returncode != 0:
        raise Exception(f"Failed to check job status: {result.stderr}")
    
    job_info = result.stdout
    job_state = None
    for line in job_info.splitlines():
        if "JobState=" in line:
            job_state = line.split("=")[1].split()[0]
            break
    return job_state

async def wait_for_job_completion(job_id):
    """Asynchronously wait for the job to complete or fail."""
    while True:
        status = parse_scontrol_job_state(job_id)
        if status in ["COMPLETED", "FAILED"]:
            return status
        await asyncio.sleep(5)  # Check every 10 seconds

def parse_output_content(output):
    try:
        # Find and evaluate the plot paths
        plot_paths_match = re.search(r"Plot Paths: (\[.*?\])", output, re.S)
        plot_paths = ast.literal_eval(plot_paths_match.group(1)) if plot_paths_match else []

        # Find and evaluate the dataset description
        dataset_desc_match = re.search(r"Dataset Description: ({.*?})", output, re.S)
        dataset_description = ast.literal_eval(dataset_desc_match.group(1)) if dataset_desc_match else {}

        return plot_paths, dataset_description

    except Exception as e:
        # Log any error that occurs during the parsing process
        print(f"Error parsing output content: {e}")
        return [], {}
    
    
counter_file = 'visit_counter.json'

def read_counter():
    try:
        with open(counter_file, "r") as file:
            data = json.load(file)
            return data['count']
    except FileNotFoundError:
        return 0

def increment_counter():
    count = read_counter()
    count += 1
    with open(counter_file, "w") as file:
        json.dump({"count": count}, file)

@app.get("/backend/visit_counter")
async def get_counter():
    increment_counter()
    return {"visit_count": read_counter()}    


app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"], 
    allow_headers=["*"],  
)