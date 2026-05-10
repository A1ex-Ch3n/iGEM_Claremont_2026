
from transformers import AutoTokenizer, EsmForMaskedLM
import torch
import csv
import pandas as pd
import numpy as np


class ESMEmbedder:

    def __init__(self, model_name: str = 'facebook/esm2_t33_650M_UR50D'):

        """
        Initializes the ESMEmbedder with a PLM model
        """

        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = EsmForMaskedLM.from_pretrained(model_name)  
        

    def _get_mean_embedding(self, sequence: str) -> np.ndarray:
        pass

    def proccess_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        pass


if __name__ == "__main__":
    pass
    
