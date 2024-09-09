from pydantic import BaseModel
from typing import Optional,List
from fastapi import Form,File,UploadFile

class UploadData(BaseModel):
    file: Optional[UploadFile] = File(None)
    species: str
    email: Optional[str] = None
    analysisName: str

class ProcessDataParams(BaseModel):
    min_genes: int = 200
    min_cells: int = 3
    mito_threshold: int = 5  
    doublet_detection: bool = True
    
class EmbeddingParams(BaseModel):
    n_neighbors: int = 15
    n_pcs: int = 18
    resolution: float = 0.5

class PcaParams(BaseModel):
    n_genes : int = 2000
    flavor : str 
    
class VisualizationInput(BaseModel):
    file_uuid: str
    dim_red: str
    gene_list: List[str]