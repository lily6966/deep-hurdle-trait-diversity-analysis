import numpy as np
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    confusion_matrix,
    precision_recall_fscore_support,
)
from sklearn.metrics import mean_squared_error, mean_absolute_error
import os
import pandas as pd
import argparse
import yaml
from utils import load_data

# ===== SETUP AND ARGUMENT PARSING =====
parser = argparse.ArgumentParser(description="Analyze hurdle model results.")
parser.add_argument(
    "--exp_dir", type=str, default="../", help="Directory containing experiment results"
)
parser.add_argument(
    "--log_dir", type=str, default="../", help="Directory containing load checkpoints testing results"
)
parser.add_argument(
    "--binary_path",
    type=str,
    default=None,
    help="Path to binary_results.npz file (will search in exp_dir if not provided)",
)
parser.add_argument(
    "--count_path",
    type=str,
    default=None,
    help="Path to count_results.npz file (will search in exp_dir if not provided)",
)
args = parser.parse_args()

# Create a results directory within the exp_dir if it doesn't exist
results_dir = os.path.join(args.exp_dir, "analysis_results")
os.makedirs(results_dir, exist_ok=True)

# ===== LOAD DATA =====
# Determine file paths
if args.binary_path:
    binary_path = args.binary_path
else:
    binary_path = os.path.join(args.exp_dir, "binary_results.npz")

if args.count_path:
    count_path = args.count_path
else:
    count_path = os.path.join(args.exp_dir, "count_results.npz")

print(f"Binary results path: {binary_path}")
print(f"Count results path: {count_path}")

# Load the data
binary_data = np.load(binary_path)["data"]
count_data = np.load(count_path)["data"]

# Split into GT and predictions
num_features = binary_data.shape[1] // 2

binary_gt = binary_data[:, :num_features]
binary_pred = binary_data[:, num_features:]

count_gt = count_data[:, :num_features]
count_pred = count_data[:, num_features:]

print(f"Binary data shape: {binary_data.shape} - {num_features} features")
print(f"Count data shape: {count_data.shape} - {num_features} features")

# ===== LOAD TRAINING DATA FOR NON-ZERO RATE CALCULATION =====
# Load hparams from the experiment directory
hparams_path = os.path.join(args.exp_dir, "hparams.yaml")
print(f"Loading hparams from {hparams_path}")
with open(hparams_path, "r") as f:
    hparams = yaml.safe_load(f)['config']

# Load training data
print("Loading training data to calculate non-zero rates...")
data = load_data(hparams)
label_names = data["label_names"]

# Calculate non-zero rate for each class
non_zero_rates = []
for i in range(num_features):
    non_zero_rate = np.mean(data['binary_label_data'][:, i] > 0)
    non_zero_rates.append(non_zero_rate)

# ===== BINARY CLASSIFICATION ANALYSIS =====
def analyze_binary_classification(
    binary_gt, binary_pred, num_features, results_dir, non_zero_rates, label_names
):
    """Analyze binary classification performance."""
    # Overall metrics
    binary_gt_flat = binary_gt.flatten()
    binary_pred_flat = np.round(binary_pred.flatten()).astype(int)

    overall_accuracy = accuracy_score(binary_gt_flat, binary_pred_flat)
    overall_f1 = f1_score(binary_gt_flat, binary_pred_flat)

    print(f"\nOverall Binary Classification:")
    print(f"Accuracy: {overall_accuracy:.4f}")
    print(f"F1 Score: {overall_f1:.4f}")

    # Confusion matrix
    cm = confusion_matrix(binary_gt_flat, binary_pred_flat)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title("Confusion Matrix")
    plt.savefig(os.path.join(results_dir, "binary_confusion_matrix.pdf"))
    plt.close()

    # Per-class analysis
    class_metrics = []

    for i in range(num_features):
        gt = binary_gt[:, i]
        pred = np.round(binary_pred[:, i]).astype(int)

        acc = accuracy_score(gt, pred)
        f1 = f1_score(gt, pred, zero_division=0)
        precision, recall, _, _ = precision_recall_fscore_support(
            gt, pred, average="binary", zero_division=0
        )
        positive_rate = np.mean(gt)

        class_metrics.append(
            {
                "Class": i,
                "Label": label_names[i],
                "Accuracy": acc,
                "F1 Score": f1,
                "Precision": precision,
                "Recall": recall,
                "Positive Rate": positive_rate,
                "Non-Zero Rate": non_zero_rates[i],
            }
        )

    df_metrics = pd.DataFrame(class_metrics)
    print("\nBinary metrics by class:")
    print(df_metrics)

    # Visualize F1 scores
    plt.figure(figsize=(6, 30))
    plt.barh(df_metrics["Label"], df_metrics["F1 Score"])
    plt.ylabel("Class")
    plt.xlabel("F1 Score")
    plt.title("F1 Score by Class")
    plt.yticks(fontsize=6)
    plt.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "binary_f1_by_class.pdf"))
    plt.close()

    # Plot non-zero rate vs F1 score
    plt.figure(figsize=(10, 8))
    plt.scatter(df_metrics["Non-Zero Rate"], df_metrics["F1 Score"], s=80)
    for i, row in df_metrics.iterrows():
        plt.annotate(
            row["Label"],
            (row["Non-Zero Rate"], row["F1 Score"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("Non-Zero Rate")
    plt.ylabel("F1 Score")
    plt.title("Non-Zero Rate vs. F1 Score")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "non_zero_rate_vs_f1.pdf"))
    plt.close()

    # Identify best and worst classes
    best_class = df_metrics.loc[df_metrics["F1 Score"].idxmax()]
    worst_class = df_metrics.loc[df_metrics["F1 Score"].idxmin()]

    print(
        f"\nBest performing class: {best_class['Label']} with F1 score of {best_class['F1 Score']:.4f}"
    )
    print(
        f"Worst performing class: {worst_class['Label']} with F1 score of {worst_class['F1 Score']:.4f}"
    )

    return overall_accuracy, overall_f1, df_metrics


# ===== COUNT REGRESSION ANALYSIS =====

def analyze_count_regression(
    count_gt, count_pred, num_features, results_dir, non_zero_rates, label_names
):
    """Analyze count regression performance."""
    # Overall RMSE
    # Convert tensors to NumPy and transform predictions from log(1 + count) scale
    if isinstance(count_pred, torch.Tensor):
        count_pred = count_pred.detach().cpu().numpy()
    if isinstance(count_gt, torch.Tensor):
        count_gt = count_gt.detach().cpu().numpy()

    #count_pred = np.expm1(count_pred)  # now in original count scale
    count_gt_flat = count_gt.flatten()
    count_pred_flat = count_pred.flatten()

    rmse = np.sqrt(mean_squared_error(count_gt_flat, count_pred_flat))
    mae = mean_absolute_error(count_gt_flat, count_pred_flat)

    print(f"\nCount Regression:")
    print(f"Overall RMSE: {rmse:.4f}")
    print(f"Overall MAE: {mae:.4f}")

    # Visualize predictions vs ground truth
    plt.figure(figsize=(10, 8))
    # Original scatter plot
    plt.subplot(1, 2, 1)
    plt.scatter(count_gt_flat, count_pred_flat, alpha=0.3)
    max_val = max(count_gt_flat.max(), count_pred_flat.max())
    plt.plot([0, max_val], [0, max_val], "r--")
    plt.xlabel("Ground Truth")
    plt.ylabel("Predictions")
    plt.title("Count Predictions vs Ground Truth")
    plt.grid(alpha=0.3)

    # Normalized scatter plot - now with per-class normalization
    plt.subplot(1, 2, 2)
    # Create normalized arrays for each class and flatten
    normalized_gt_list = []
    normalized_pred_list = []

    for i in range(num_features):
        gt = count_gt[:, i]
        pred = count_pred[:, i]

        max_val = max(gt.max(), 1)  # Use at least 1 to avoid division by zero

        # Normalize by class max
        normalized_gt_list.append(gt / max_val)
        normalized_pred_list.append(pred / max_val)

    # Flatten for scatter plot
    normalized_gt = np.concatenate(normalized_gt_list)
    normalized_pred = np.concatenate(normalized_pred_list)

    plt.scatter(normalized_gt, normalized_pred, alpha=0.3)
    plt.plot([0, 1], [0, 1], "r--")
    plt.xlabel("Normalized Ground Truth")
    plt.ylabel("Normalized Predictions")
    plt.title("Per-Class Normalized Count Predictions (0-1)")
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "count_predictions.pdf"))
    plt.close()

    # Per-class RMSE analysis
    class_rmse = []

    for i in range(num_features):
        gt = count_gt[:, i]
        pred = count_pred[:, i]

        rmse_val = np.sqrt(mean_squared_error(gt, pred))
        mae_val = mean_absolute_error(gt, pred)
        avg_count = np.mean(gt)
        max_count = np.max(gt)

        # Calculate normalized RMSE (as percentage of average)
        nrmse_avg = rmse_val / avg_count if avg_count > 0 else 0

        # Calculate normalized RMSE (as percentage of max)
        nrmse_max = rmse_val / max_count if max_count > 0 else 0

        # Normalize by class max
        class_max = max(gt.max(), 1)  # Use at least 1 to avoid division by zero
        normalized_gt = gt / class_max
        normalized_pred = pred / class_max
        normalized_rmse = np.sqrt(mean_squared_error(normalized_gt, normalized_pred))

        class_rmse.append(
            {
                "Class": i,
                "Label": label_names[i],
                "RMSE": rmse_val,
                "MAE": mae_val,
                "NRMSE (by avg)": nrmse_avg,
                "NRMSE (by max)": nrmse_max,
                "Normalized RMSE (0-1)": normalized_rmse,
                "Avg Count": avg_count,
                "Max Count": max_count,
                "Non-Zero Rate": non_zero_rates[i],
            }
        )

    df_rmse = pd.DataFrame(class_rmse)
    print("\nRMSE by class:")
    print(df_rmse)

    # Visualize RMSE by class
    plt.figure(figsize=(14, 8))

    # Original RMSE
    plt.subplot(2, 1, 1)
    plt.bar(df_rmse["Label"], df_rmse["RMSE"])
    plt.xlabel("Class")
    plt.ylabel("RMSE")
    plt.title("RMSE by Class")
    plt.xticks(rotation=90)
    plt.grid(axis="y", alpha=0.3)

    # Normalized RMSE (0-1 scale)
    plt.subplot(2, 1, 2)
    plt.bar(df_rmse["Label"], df_rmse["Normalized RMSE (0-1)"])
    plt.xlabel("Class")
    plt.ylabel("Normalized RMSE")
    plt.title("Normalized RMSE by Class (0-1 scale)")
    plt.xticks(rotation=90)
    plt.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "count_rmse_by_class.pdf"))
    plt.close()




    # Plot non-zero rate vs RMSE
    plt.figure(figsize=(12, 10))

    # Original RMSE vs non-zero rate
    plt.subplot(2, 1, 1)
    plt.scatter(df_rmse["Non-Zero Rate"], df_rmse["RMSE"], s=80)
    for i, row in df_rmse.iterrows():
        plt.annotate(
            row["Label"],
            (row["Non-Zero Rate"], row["RMSE"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("Non-Zero Rate")
    plt.ylabel("RMSE")
    plt.title("Non-Zero Rate vs. RMSE")
    plt.grid(alpha=0.3)

    # Normalized RMSE vs non-zero rate
    plt.subplot(2, 1, 2)
    plt.scatter(df_rmse["Non-Zero Rate"], df_rmse["Normalized RMSE (0-1)"], s=80)
    for i, row in df_rmse.iterrows():
        plt.annotate(
            row["Label"],
            (row["Non-Zero Rate"], row["Normalized RMSE (0-1)"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("Non-Zero Rate")
    plt.ylabel("Normalized RMSE (0-1)")
    plt.title("Non-Zero Rate vs. Normalized RMSE")
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "non_zero_rate_vs_rmse.pdf"))
    plt.close()

    # Identify best and worst classes
    best_class = df_rmse.loc[df_rmse["RMSE"].idxmin()]
    worst_class = df_rmse.loc[df_rmse["RMSE"].idxmax()]

    print(
        f"\nBest performing class: {best_class['Label']} with RMSE of {best_class['RMSE']:.4f}"
    )
    print(
        f"Worst performing class: {worst_class['Label']} with RMSE of {worst_class['RMSE']:.4f}"
    )

    return rmse, mae, df_rmse


# ===== HURDLE MODEL ANALYSIS =====
def analyze_hurdle_model(
    binary_gt,
    binary_pred,
    count_gt,
    count_pred,
    df_metrics,
    df_rmse,
    num_features,
    results_dir,
    label_names,
):
    """Analyze combined hurdle model performance."""
    # Combine binary and count metrics
    combined = pd.merge(df_metrics, df_rmse, on=["Class", "Label", "Non-Zero Rate"])
    print("\nCombined metrics:")
    print(
        combined[
            ["Label", "F1 Score", "RMSE", "Positive Rate", "Avg Count", "Non-Zero Rate"]
        ]
    )

    # Analyze relationship between binary and count performance
    plt.figure(figsize=(12, 6))

    # Original plot
    plt.subplot(1, 2, 1)
    plt.scatter(combined["F1 Score"], combined["RMSE"], s=80)
    # Label points with class labels
    for i, row in combined.iterrows():
        plt.annotate(
            row["Label"],
            (row["F1 Score"], row["RMSE"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("F1 Score (Binary)")
    plt.ylabel("RMSE (Count)")
    plt.title("Binary vs Count Performance")
    plt.grid(alpha=0.3)

    # Normalized plot
    plt.subplot(1, 2, 2)
    plt.scatter(combined["F1 Score"], combined["Normalized RMSE (0-1)"], s=80)
    # Label points with class labels
    for i, row in combined.iterrows():
        plt.annotate(
            row["Label"],
            (row["F1 Score"], row["Normalized RMSE (0-1)"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("F1 Score (Binary)")
    plt.ylabel("Normalized RMSE (0-1)")
    plt.title("Binary vs Normalized Count Performance")
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "binary_vs_count_performance.pdf"))
    plt.close()

    # Plot with non-zero rate as color
    plt.figure(figsize=(12, 6))

    # Original plot with non-zero rate as color
    plt.subplot(1, 2, 1)
    scatter = plt.scatter(
        combined["F1 Score"],
        combined["RMSE"],
        c=combined["Non-Zero Rate"],
        cmap="viridis",
        s=100,
        alpha=0.8,
    )
    # Label points with class labels
    for i, row in combined.iterrows():
        plt.annotate(
            row["Label"],
            (row["F1 Score"], row["RMSE"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("F1 Score (Binary)")
    plt.ylabel("RMSE (Count)")
    plt.title("Binary vs Count Performance (colored by Non-Zero Rate)")
    plt.colorbar(scatter, label="Non-Zero Rate")
    plt.grid(alpha=0.3)

    # Normalized plot with non-zero rate as color
    plt.subplot(1, 2, 2)
    scatter = plt.scatter(
        combined["F1 Score"],
        combined["Normalized RMSE (0-1)"],
        c=combined["Non-Zero Rate"],
        cmap="viridis",
        s=100,
        alpha=0.8,
    )
    # Label points with class labels
    for i, row in combined.iterrows():
        plt.annotate(
            row["Label"],
            (row["F1 Score"], row["Normalized RMSE (0-1)"]),
            xytext=(5, 0),
            textcoords="offset points",
        )
    plt.xlabel("F1 Score (Binary)")
    plt.ylabel("Normalized RMSE (0-1)")
    plt.title("Binary vs Normalized Count Performance (colored by Non-Zero Rate)")
    plt.colorbar(scatter, label="Non-Zero Rate")
    plt.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "binary_vs_count_performance_with_nzr.pdf"))
    plt.close()

    # Analyze impact of binary correctness on count prediction
    all_binary_gt = binary_gt.flatten()
    all_binary_pred = np.round(binary_pred.flatten()).astype(int)
    all_count_gt = count_gt.flatten()
    all_count_pred = count_pred.flatten()

    # Keep track of which feature each value belongs to
    feature_indices = np.repeat(np.arange(num_features), count_gt.shape[0])

    # Identify where binary prediction was correct/incorrect
    binary_correct = all_binary_gt == all_binary_pred
    binary_incorrect = ~binary_correct

    # Calculate RMSE for when binary prediction was correct vs incorrect
    if np.sum(binary_correct) > 0:
        rmse_when_correct = np.sqrt(
            mean_squared_error(
                all_count_gt[binary_correct], all_count_pred[binary_correct]
            )
        )
        print(f"\nRMSE when binary prediction is correct: {rmse_when_correct:.4f}")

        # Normalized RMSE when correct - using per-class normalization
        correct_gt = []
        correct_pred = []

        # Get the feature indices for correct predictions
        correct_features = feature_indices[binary_correct]

        # Group by feature and normalize
        for i in range(num_features):
            feature_mask = correct_features == i
            if np.sum(feature_mask) > 0:
                gt_values = all_count_gt[binary_correct][feature_mask]
                pred_values = all_count_pred[binary_correct][feature_mask]

                # Get max for this feature from the original data
                feature_max = max(count_gt[:, i].max(), 1)

                # Normalize
                correct_gt.extend(gt_values / feature_max)
                correct_pred.extend(pred_values / feature_max)

        if correct_gt:
            norm_rmse_correct = np.sqrt(
                mean_squared_error(np.array(correct_gt), np.array(correct_pred))
            )
            print(
                f"Per-class normalized RMSE when binary prediction is correct: {norm_rmse_correct:.4f}"
            )

    if np.sum(binary_incorrect) > 0:
        rmse_when_incorrect = np.sqrt(
            mean_squared_error(
                all_count_gt[binary_incorrect], all_count_pred[binary_incorrect]
            )
        )
        print(f"RMSE when binary prediction is incorrect: {rmse_when_incorrect:.4f}")

        # Normalized RMSE when incorrect - using per-class normalization
        incorrect_gt = []
        incorrect_pred = []

        # Get the feature indices for incorrect predictions
        incorrect_features = feature_indices[binary_incorrect]

        # Group by feature and normalize
        for i in range(num_features):
            feature_mask = incorrect_features == i
            if np.sum(feature_mask) > 0:
                gt_values = all_count_gt[binary_incorrect][feature_mask]
                pred_values = all_count_pred[binary_incorrect][feature_mask]

                # Get max for this feature from the original data
                feature_max = max(count_gt[:, i].max(), 1)

                # Normalize
                incorrect_gt.extend(gt_values / feature_max)
                incorrect_pred.extend(pred_values / feature_max)

        if incorrect_gt:
            norm_rmse_incorrect = np.sqrt(
                mean_squared_error(np.array(incorrect_gt), np.array(incorrect_pred))
            )
            print(
                f"Per-class normalized RMSE when binary prediction is incorrect: {norm_rmse_incorrect:.4f}"
            )

    return combined


# ===== RUN ANALYSES =====
# Binary Classification Analysis
overall_accuracy, overall_f1, df_metrics = analyze_binary_classification(
    binary_gt, binary_pred, num_features, results_dir, non_zero_rates, label_names
)

# Count Regression Analysis
rmse, mae, df_rmse = analyze_count_regression(
    count_gt, count_pred, num_features, results_dir, non_zero_rates, label_names
)

# Hurdle Model Analysis
combined = analyze_hurdle_model(
    binary_gt,
    binary_pred,
    count_gt,
    count_pred,
    df_metrics,
    df_rmse,
    num_features,
    results_dir,
    label_names,
)

# ===== SAVE RESULTS AND SUMMARY =====
# Save metrics to CSV files
df_metrics.to_csv(os.path.join(results_dir, "binary_metrics.csv"), index=False)
df_rmse.to_csv(os.path.join(results_dir, "count_metrics.csv"), index=False)
combined.to_csv(os.path.join(results_dir, "combined_metrics.csv"), index=False)

# Save non-zero rates specifically
non_zero_df = pd.DataFrame(
    {
        "Class": range(num_features),
        "Label": label_names,
        "Non-Zero Rate": non_zero_rates,
    }
)
non_zero_df.to_csv(os.path.join(results_dir, "non_zero_rates.csv"), index=False)

# Summary
print("\nSummary and Conclusions:")
print("1. Binary Classification:")
print(f"   - Overall accuracy: {overall_accuracy:.4f}, F1 score: {overall_f1:.4f}")
best_class_idx = df_metrics["F1 Score"].idxmax()
worst_class_idx = df_metrics["F1 Score"].idxmin()
print(
    f"   - Best class: {df_metrics.loc[best_class_idx]['Label']} (F1={df_metrics['F1 Score'].max():.4f})"
)
print(
    f"   - Worst class: {df_metrics.loc[worst_class_idx]['Label']} (F1={df_metrics['F1 Score'].min():.4f})"
)

print("\n2. Count Prediction:")
print(f"   - Overall RMSE: {rmse:.4f}")
best_class_idx = df_rmse["RMSE"].idxmin()
worst_class_idx = df_rmse["RMSE"].idxmax()
print(
    f"   - Best class: {df_rmse.loc[best_class_idx]['Label']} (RMSE={df_rmse['RMSE'].min():.4f})"
)
print(
    f"   - Worst class: {df_rmse.loc[worst_class_idx]['Label']} (RMSE={df_rmse['RMSE'].max():.4f})"
)

print("\n3. Per-Class Normalized Count Prediction (0-1 scale):")
best_class_idx = df_rmse["Normalized RMSE (0-1)"].idxmin()
worst_class_idx = df_rmse["Normalized RMSE (0-1)"].idxmax()
print(
    f"   - Best class: {df_rmse.loc[best_class_idx]['Label']} (Normalized RMSE={df_rmse['Normalized RMSE (0-1)'].min():.4f})"
)
print(
    f"   - Worst class: {df_rmse.loc[worst_class_idx]['Label']} (Normalized RMSE={df_rmse['Normalized RMSE (0-1)'].max():.4f})"
)

# Add a note about the normalization method
print(
    "\nNote: Count values were normalized by dividing by each class's maximum value to create"
)
print(
    "a 0-1 scale for each class independently, allowing for fair comparisons between features with vastly different magnitudes."
)

print(f"\nAnalysis results saved to: {results_dir}")
