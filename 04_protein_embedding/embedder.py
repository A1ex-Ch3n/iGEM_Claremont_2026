
from transformers import AutoTokenizer, EsmModel
import torch
import csv
import pandas as pd
import numpy as np
from tqdm import tqdm


class ESMEmbedder:

    def __init__(self, model_name: str = 'facebook/esm2_t33_650M_UR50D'):

        """
        Initializes the ESMEmbedder with a PLM model
        """

        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = EsmModel.from_pretrained(model_name)  
        


    def _get_mean_embedding(self, sequences: list[str]) -> np.ndarray:

        """
        Generates the mean embeddings of a protein sequence
        """

        inputs = self.tokenizer(sequences, return_tensors="pt", padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs, output_hidden_states=True)
            last_hidden_state = outputs.hidden_states[-1]
            attention_mask = inputs["attention_mask"]


            embedding = (last_hidden_state * attention_mask.unsqueeze(-1)).sum(1) / attention_mask.sum(1).unsqueeze(-1)
            return embedding.cpu().numpy()

    def proccess_dataframe(self, df: pd.DataFrame, seq_col: str, batch_size: int = 32) -> pd.DataFrame:


        """
        Batch process protein sequences for embedding
        """

        sequences = df[seq_col].tolist()
        all_embeddings = []

        print(f"Processing {len(sequences)} sequences in batches of {batch_size}...")

        for i in tqdm(range(0, len(sequences), batch_size)):
        
            batch_seqs = sequences[i : i + batch_size]

            batch_embeddings = self._get_mean_embedding(batch_seqs)

            all_embeddings.extend(batch_embeddings)

        df["embeddings"] = all_embeddings

        return df


if __name__ == "__main__":
    pass
    
