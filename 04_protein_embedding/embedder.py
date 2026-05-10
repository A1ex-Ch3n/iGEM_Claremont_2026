
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

        """
        Generates the mean embeddings of a protein sequence
        """


        inputs = self.tokenizer(sequence, return_tensors="pt", padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs, output_hidden_states=True)
            last_hidden_state = outputs.hidden_states[-1]
            attention_mask = inputs["attention_mask"]
            embedding = (last_hidden_state * attention_mask.unsqueeze(-1)).sum(1) / attention_mask.sum(1).unsqueeze(-1)
            return embedding.squeeze(0).numpy()

    def proccess_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        pass


if __name__ == "__main__":
    pass
    
