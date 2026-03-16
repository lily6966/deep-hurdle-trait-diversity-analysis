# Deep Hurdle Population Model

A PyTorch implementation of a deep hurdle model for population modeling.

## Overview

This model implements a two-part hurdle model with:
1. A classifier that predicts presence/absence (binary outcome)
2. A regressor that predicts counts (continuous outcome)

The model is implemented using PyTorch Lightning for efficient training and evaluation.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/Deep-hurdle-population-model.git
cd Deep-hurdle-population-model
```

2. Create a virtual environment and install dependencies:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Data Preparation

The model expects data in a specific format. Your data should be stored as an `.npz` file with the following structure:
- Features array
- Count labels array
- Binary presence/absence labels array

Update the `data_path` in your configuration file (`config.yaml`) to point to your data file.

## Training

### Basic Training

To train the model using the default configuration:

```bash
python train.py
```

This will:
1. Load configuration from `config.yaml`
2. Train the classifier for the first half of epochs
3. Train the regressor for the second half of epochs
4. Evaluate on the test set
5. Run analysis on the results

### Custom Configuration

To use a custom configuration file:

```bash
python train.py --config path/to/your/config.yaml
```

### Named Experiments

To name your experiment for easier tracking:

```bash
python train.py --exp experiment_name
```

Experiment results will be saved in `experiments/experiment_name/version_X/`.

### Configuration Options

Edit `config.yaml` to customize the model settings:

```yaml
# Model configuration
model:
  input_dim: 33             # Input feature dimension
  label_dim: 395            # Output label dimension
  latent_dim: 64            # Latent space dimension
  emb_size: 2048            # Embedding size
  scale_coeff: 1            # Scale coefficient
  reg: "gmvae"              # Regularization type
  test_sample: false        # Whether to sample during testing
  keep_prob: 0.5            # Dropout keep probability

# Training configuration
training:
  max_epoch: 50             # Maximum number of epochs per stage
  batch_size: 4096          # Batch size
  lr_classifier: 5e-3       # Learning rate for classifier
  lr_regressor: 2e-4        # Learning rate for regressor
  lr_decay: 0.5             # Learning rate decay factor
  lr_decay_times: 4         # Number of times to decay learning rate

# Data configuration
data:
  data_path: './DATA/ebird_entire/ebird_data.npz'  # Path to data file
  train_size: 0.7           # Proportion of data for training
  val_size: 0.15            # Proportion of data for validation
  test_size: 0.15           # Proportion of data for testing
  random_seed: 42           # Random seed for reproducibility
```

## Testing Trained Models

To evaluate a trained model:

```bash
python train.py --ckpt path/to/checkpoint.ckpt
```

This will:
1. Load the model from the checkpoint
2. Run evaluation on the test set
3. Generate analysis reports

## Analysis

After training or testing, the model automatically runs the analysis script, which:
1. Calculates classification metrics (accuracy, precision, recall, F1)
2. Calculates regression metrics (RMSE, MAE)
3. Generates visualizations
4. Saves results to the experiment directory

You can also run the analysis separately:

```bash
python analyze_results.py --exp_dir experiments/experiment_name/version_X
```

## Output

Training output is organized as follows:
- `experiments/experiment_name/version_X/checkpoints/`: Model checkpoints
- `experiments/experiment_name/version_X/analysis_log.txt`: Analysis results
- TensorBoard logs for tracking metrics

To view training metrics, run:
```bash
tensorboard --logdir experiments/experiment_name
``` 