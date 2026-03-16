import torch
import numpy as np
import yaml
from sklearn.model_selection import train_test_split


class Dataset(torch.utils.data.Dataset):
    """Legacy dataset class (kept for backward compatibility)"""

    def __init__(self, features, orig_labels, log_count_labels, binary_labels):
        
        self.features = features
        self.orig_labels = orig_labels
        self.log_count_labels = log_count_labels
        self.binary_labels = binary_labels

    def __len__(self):
        return len(self.features)

    def __getitem__(self, index):
        #loc = np.array(self.loc[index], dtype=np.float32)
        X = self.features[index]
        count = self.orig_labels[index]
        log_count = self.log_count_labels[index]
        binary = self.binary_labels[index]
        return (X, count, log_count, binary)


def load_config(config_path="config.yaml"):
    """Load configuration from a YAML file"""
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)
    return config


def load_esrd():
    """Load ESRD data from a specified path"""
    esrd_data = np.load('./DATA/ebird_entire/ebird_data_esrd0718.npz', allow_pickle=True)
    
    # Extract features and labels
    feat_data = esrd_data["feat"]
    #loc = esrd_data["loc"]
    orig_label_data = esrd_data["count_label"]
    
    # Generate binary labels from count labels
    binary_label_data = (orig_label_data > 0) 
    # Log-transform labels
    count_label_data = np.log1p(orig_label_data)

    # Store data in dictionary
    data = {
        #"loc": loc,
        "feature_names": esrd_data["feature_names"],
        "label_names": esrd_data["label_names"],
        "feat_data": feat_data.astype(np.float32),
        "orig_count_label": orig_label_data.astype(np.float32),
        "count_label_data": count_label_data.astype(np.float32),
        "binary_label_data": binary_label_data.astype(np.float32),

    }

    return data

def load_esrd_low():
    """Load ESRD data from a specified path"""
    esrd_data = np.load('./DATA/ebird_entire/ebird_data_esrd_low.npz', allow_pickle=True)
    
    # Extract features and labels
    feat_data = esrd_data["feat"]
    #loc = esrd_data["loc"]
    orig_label_data = esrd_data["count_label"]
    
    # Generate binary labels from count labels
    binary_label_data = (orig_label_data > 0) 
    # Log-transform labels
    count_label_data = np.log1p(orig_label_data)

    # Store data in dictionary
    data = {
        #"loc": loc,
        "feature_names": esrd_data["feature_names"],
        "label_names": esrd_data["label_names"],
        "feat_data": feat_data.astype(np.float32),
        "orig_count_label": orig_label_data.astype(np.float32),
        "count_label_data": count_label_data.astype(np.float32),
        "binary_label_data": binary_label_data.astype(np.float32),

    }

    return data


def load_esrd_inthigh():
    """Load ESRD data from a specified path"""
    esrd_data = np.load('./DATA/ebird_entire/ebird_data_esrd_inthigh.npz', allow_pickle=True)
    
    # Extract features and labels
    feat_data = esrd_data["feat"]
    #loc = esrd_data["loc"]
    orig_label_data = esrd_data["count_label"]
    
    # Generate binary labels from count labels
    binary_label_data = (orig_label_data > 0) 
    # Log-transform labels
    count_label_data = np.log1p(orig_label_data)

    # Store data in dictionary
    data = {
        #"loc": loc,
        "feature_names": esrd_data["feature_names"],
        "label_names": esrd_data["label_names"],
        "feat_data": feat_data.astype(np.float32),
        "orig_count_label": orig_label_data.astype(np.float32),
        "count_label_data": count_label_data.astype(np.float32),
        "binary_label_data": binary_label_data.astype(np.float32),

    }

    return data

def load_data(config):
    """Load and preprocess data from the NPZ file"""

    # Load data from NPZ file
    data_file = np.load(config["data"]["data_path"], allow_pickle=True)

    # Extract features and labels
    feat_data = data_file["feat"]
    #loc = data_file["loc"]
    count_label_data = data_file["count_label"]
    # Generate binary labels from count labels
    binary_label_data = (count_label_data > 0) 
    # Log-transform labels
    log_count_label = np.log1p(count_label_data)
    
    


    
    feat_data = np.array(feat_data)
    print("feat_data type:", type(feat_data))
    print("feat_data shape:", getattr(feat_data, 'shape', 'no shape'))

    # 🔍 Add this line to debug 'np'
    print("np is actually:", np, "| type:", type(np))


    # Normalize features
    feat_data = (feat_data - np.mean(feat_data, axis=0)) / (
        np.std(feat_data, axis=0) + 1e-8
    )
    

    # Split data into train/val/test
    indices = np.arange(len(feat_data))

    # First split into train and temp (val+test)
    train_idx, temp_idx = train_test_split(
        indices,
        train_size=config["data"]["train_size"],
        random_state=config["data"]["random_seed"],
    )

    # Then split temp into val and test
    relative_val_size = config["data"]["val_size"] / (1 - config["data"]["train_size"])
    val_idx, test_idx = train_test_split(
        temp_idx,
        train_size=relative_val_size,
        random_state=config["data"]["random_seed"],
    )

    # Store data in dictionary
    data = {
        #"loc": loc,
        "feature_names": data_file["feature_names"],
        "label_names": data_file["label_names"],
        "feat_data": feat_data.astype(np.float32),
        "orig_count_label": count_label_data.astype(np.float32),
        "count_label_data": log_count_label.astype(np.float32),
        "binary_label_data": binary_label_data.astype(np.float32),
        "train_idx": train_idx,
        "val_idx": val_idx,
        "test_idx": test_idx,
    }

    return data


def get_dataloaders(data, config, split="train"):
    """Create data loaders for training, validation, or testing"""
    
    batch_size = config["training"]["batch_size"]
    shuffle = split == "train"
    idx_key = f"{split}_idx"

    dataset = Dataset(
        #data["loc"][data[idx_key]],
        data["feat_data"][data[idx_key]],
        data["orig_count_label"][data[idx_key]],
        data["count_label_data"][data[idx_key]],
        data["binary_label_data"][data[idx_key]],
    )

    return torch.utils.data.DataLoader(
        dataset, batch_size=batch_size, shuffle=shuffle,
        num_workers=0,
    )


def get_dataloaders_esrd(data, config):
    """Create data loaders for ESRD data"""
    
    batch_size = config["training"]["batch_size"]
    shuffle = False  # ESRD data is not shuffled

    dataset = Dataset(
        #data["loc"],
        data["feat_data"],
        data["orig_count_label"],
        data["count_label_data"],
        data["binary_label_data"],
    )

    return torch.utils.data.DataLoader(
        dataset, batch_size=batch_size, shuffle=shuffle,
        num_workers=0,
    )