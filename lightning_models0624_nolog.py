import torch
import torch.nn as nn
import pytorch_lightning as pl
from torch.optim import Adam
from torch.optim.lr_scheduler import StepLR
from model import VAE, compute_loss
import torch.nn.functional as F
import os
import numpy as np
import torchmetrics
import pandas as pd

class DeepHurdleModule(pl.LightningModule):
    def __init__(self, config):
        super().__init__()
        self.save_hyperparameters()
        self.config = config

        # Convert config to args for backward compatibility
        self.classifier_args = self._config_to_args(config, mode="classification")
        self.regressor_args = self._config_to_args(config, mode="regression")

        # Initialize models
        self.classifier = VAE(self.classifier_args)
        self.regressor = VAE(self.regressor_args)

        # Define loss functions
        self.classifier_criterion = nn.BCELoss()
        self.regressor_criterion = nn.PoissonNLLLoss(log_input=False, reduction="none")

        self.classifier_trained = False
        self.automatic_optimization = False

        # Lists to store esrd predictions and ground truth
        self.test_binary_gt = []
        self.test_count_gt = []
        self.test_binary_pred = []
        self.test_count_pred = []
        self.test_loc = []  # Store locations for later use
        self.test_label_regressor = []
        self.test_label_classifier = []
        self.test_feature_classifier = []
        self.test_feature_regressor = []

        # Define metrics using torchmetrics
        self.train_accuracy = torchmetrics.Accuracy(task="binary")
        self.val_accuracy = torchmetrics.Accuracy(task="binary")
        self.test_accuracy = torchmetrics.Accuracy(task="binary")

        self.train_mse = torchmetrics.MeanSquaredError()
        self.val_mse = torchmetrics.MeanSquaredError()
        self.test_mse = torchmetrics.MeanSquaredError()


    def _config_to_args(self, config, mode="classification"):
        """Convert config dict to args object for backward compatibility"""

        class Args:
            pass

        args = Args()

        # Model parameters
        args.input_dim = config["model"]["input_dim"]
        args.label_dim = config["model"]["label_dim"]
        args.latent_dim = config["model"]["latent_dim"]
        args.emb_size = config["model"]["emb_size"]
        args.scale_coeff = config["model"]["scale_coeff"]
        args.reg = config["model"]["reg"]
        args.test_sample = config["model"]["test_sample"]
        args.keep_prob = config["model"]["keep_prob"]

        # Training parameters
        args.max_epoch = config["training"]["max_epoch"]
        args.batch_size = config["training"]["batch_size"]

        if mode == "classification":
            args.lr = config["training"]["lr_classifier"]
        else:
            args.lr = config["training"]["lr_regressor"]

        args.lr_decay = config["training"]["lr_decay"]
        args.lr_decay_times = config["training"]["lr_decay_times"]
        args.mode = mode

        return args

    def on_train_start(self):
        """Reset classifier_trained flag at the start of training"""
        self.classifier_trained = False

    def on_train_epoch_start(self):
        """Set classifier_trained flag if we're past halfway point in training"""
        if self.trainer.current_epoch > self.classifier_args.max_epoch:
            self.classifier_trained = True
            self.classifier.eval()
            self.regressor.train()
        else:
            self.classifier.train()
            self.regressor.eval()

    def on_train_epoch_end(self):
        """Update learning rate after each epoch"""
        if not self.classifier_trained:
            scheduler = self.lr_schedulers()[0]
        else:
            scheduler = self.lr_schedulers()[1]
        scheduler.step()

        # Log accumulated metrics at the end of the epoch
        if not self.classifier_trained:
            self.log("train/accuracy", self.train_accuracy.compute())
        else:
            self.log("train/rmse", torch.sqrt(self.train_mse.compute()))

        # Reset metrics for next epoch
        self.train_accuracy.reset()
        self.train_mse.reset()

    def training_step(self, batch, batch_idx):
        """
        Training step for either classifier or regressor

        Args:
            batch: Batch data containing features, count labels, and binary labels
            batch_idx: Batch index
            optimizer_idx: 0 for classifier, 1 for regressor

        Returns:
            Loss dictionary
        """
        features, orig_labels, count_labels, binary_labels = batch
        
        # get device from inputs (usually your model/device)
        if not self.classifier_trained:
            optimizer = self.optimizers()[0]
        else:
            optimizer = self.optimizers()[1]

        optimizer.zero_grad()

        # if epoch > half of max_epoch, switch to regressor training
        if not self.classifier_trained:
            # Train classifier
            output = self.classifier(features, orig_labels, binary_labels)
            total_loss, nll_loss, nll_loss_x, kl_loss, cpc_loss, pred_e, pred_x = (
                compute_loss(
                    binary_labels,
                    orig_labels,
                    output,
                    self.classifier_criterion,
                    mode="classification",
                )
            )

            # Compute column-wise mean
            #mean_per_col = pred_x.mean(dim=0, keepdim=True)  # shape [1, label_dim]

            # Binarize: 1 if greater than column mean, else 0
            #pred_binary_rounded = (pred_x > mean_per_col).float()

            # Calculate accuracy using torchmetrics
            mean_per_col = pred_x.mean(dim=0)
            #quantile_25 = pred_x_classifier.quantile(0.25, dim=0)             # shape [S]

            # Masks
            high_mean_mask = (mean_per_col >= 0.2)                 # Rule 1
            mid_mean_mask = (mean_per_col < 0.2) & (mean_per_col > 0.01)  # Rule 2
            low_mean_mask = mean_per_col <= 0.01                   # Rule 3

            # Expand masks to match pred_x shape [N, S]
            high_mask_exp = high_mean_mask.unsqueeze(0).expand_as(pred_x)
            mid_mask_exp = mid_mean_mask.unsqueeze(0).expand_as(pred_x)
            low_mask_exp = low_mean_mask.unsqueeze(0).expand_as(pred_x)

            # Thresholds
            bin_high = (pred_x > mean_per_col).float()                      # Rule 1 threshold
            bin_mid = (pred_x > mean_per_col).float()               # Rule 2 threshold
            bin_low = (pred_x > 0.05).float()                     # Rule 3 all zeros

            # Apply conditional logic
            pred_binary_rounded = torch.where(
                high_mask_exp, bin_high,
                torch.where(
                    mid_mask_exp, bin_mid,
                    bin_low
                )
            )
            self.train_accuracy(pred_binary_rounded, binary_labels)

            self.log("train/classifier_loss", total_loss, prog_bar=True)
        else:
            # Train regressor with classifier predictions
            with torch.no_grad():
                classifier_output = self.classifier(
                    features, orig_labels, binary_labels
                )
                pred_binary = torch.sigmoid(classifier_output["feat_out"])

            # Forward pass through regressor
            output = self.regressor(features, orig_labels, binary_labels)
            total_loss, nll_loss, nll_loss_x, kl_loss, cpc_loss, pred_e, pred_x = (
                compute_loss(
                    binary_labels,
                    orig_labels,
                    output,
                    self.regressor_criterion,
                    mode="regression",
                    pred_binary=pred_binary,
                )
            )

            # Calculate RMSE metrics
            # If pred_binary < 0.01, set pred_x from regressor as 0, otherwise keep original
            adjusted_pred_x = torch.where(
                pred_binary < 0.01, torch.zeros_like(pred_x), pred_x
            )
            # Inverse transform predictions and labels
            pred_x_real = adjusted_pred_x


            # Round pred_x to integers
            rounded_pred_x = torch.round(pred_x_real)
            
            #device = torch.device("mps" if torch.has_mps else "cpu")
            # Move both to CPU before passing to metric
            rounded_pred_x = rounded_pred_x.cpu()
            orig_labels = orig_labels.cpu()
            
            self.train_mse(rounded_pred_x, orig_labels)

            self.log("train/regressor_loss", total_loss, prog_bar=True)

        self.manual_backward(total_loss)
        optimizer.step()

    def validation_step(self, batch, batch_idx):
        """
        Validation step for either classifier or regressor
        """
        features, orig_labels, count_labels, binary_labels = batch
        
        # get device from inputs (usually your model/device)

        classifier_output = self.classifier(features, orig_labels, binary_labels)
        regressor_output = self.regressor(features, orig_labels, binary_labels)

        classifier_loss, _, _, _, _, _, pred_x_classifier = compute_loss(
            binary_labels,
            orig_labels,
            classifier_output,
            self.classifier_criterion,
            mode="classification",
        )

        pred_binary = torch.sigmoid(classifier_output["feat_out"])

        regressor_loss, _, _, _, _, _, pred_x_regressor = compute_loss(
            binary_labels,
            orig_labels,
            regressor_output,
            self.regressor_criterion,
            mode="regression",
            pred_binary=pred_binary,
        )

        # Calculate accuracy using torchmetrics
        # Calculate accuracy using torchmetrics
        mean_per_col = pred_x_classifier.mean(dim=0)
        #quantile_25 = pred_x_classifier.quantile(0.25, dim=0)             # shape [S]

        # Masks
        high_mean_mask = (mean_per_col >= 0.2)                 # Rule 1
        mid_mean_mask = (mean_per_col < 0.2) & (mean_per_col > 0.01)  # Rule 2
        low_mean_mask = mean_per_col <= 0.01                   # Rule 3

        # Expand masks to match pred_x shape [N, S]
        high_mask_exp = high_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)
        mid_mask_exp = mid_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)
        low_mask_exp = low_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)

        # Thresholds
        bin_high = (pred_x_classifier > mean_per_col).float()                      # Rule 1 threshold
        bin_mid = (pred_x_classifier > mean_per_col).float()               # Rule 2 threshold
        bin_low = (pred_x_classifier > 0.05).float()                     # Rule 3 all zeros

        # Apply conditional logic
        pred_binary_rounded = torch.where(
            high_mask_exp, bin_high,
            torch.where(
                mid_mask_exp, bin_mid,
                bin_low
            )
        )
        self.val_accuracy(pred_binary_rounded, binary_labels)

        # Calculate RMSE using torchmetrics
        # If pred_binary < 0.01, set pred_x from regressor as 0, otherwise keep original
        adjusted_pred_x = torch.where(
            pred_binary < 0.01, torch.zeros_like(pred_x_regressor), pred_x_regressor
        )
      
        # Inverse transform predictions and labels
        pred_x_real = adjusted_pred_x


        # Round pred_x to integers
        rounded_pred_x = torch.round(pred_x_real)
        #device = torch.device("mps" if torch.has_mps else "cpu")
        #orig_labels = orig_labels.to(device)
        #self.val_mse = self.val_mse.to(self.device)  # <--- move metric to device here
        # Move labels to device
        #orig_labels = orig_labels.to(device)
        
        # Move both to CPU before passing to metric
        rounded_pred_x = rounded_pred_x.cpu()
        orig_labels = orig_labels.cpu()
        self.val_mse(rounded_pred_x, orig_labels)
       
        # Log per-batch loss metrics
        self.log("val/classifier_loss", classifier_loss)
        self.log("val/regressor_loss", regressor_loss)

    def on_validation_epoch_end(self):
        """Log accumulated validation metrics at the end of the epoch"""
        self.log("val/accuracy", self.val_accuracy.compute())
        self.log("val/rmse", torch.sqrt(self.val_mse.compute()))

        # Reset metrics for next epoch
        self.val_accuracy.reset()
        self.val_mse.reset()

    def test_step(self, batch, batch_idx):
        """
        Test step to collect predictions and ground truth
        """
        features, orig_labels, count_labels, binary_labels = batch

        classifier_output = self.classifier(features, orig_labels, binary_labels)
        regressor_output = self.regressor(features, orig_labels, binary_labels)

        classifier_loss, _, _, _, _, _, pred_x_classifier = compute_loss(
            binary_labels,
            orig_labels,
            classifier_output,
            self.classifier_criterion,
            mode="classification",
        )

        pred_binary = torch.sigmoid(classifier_output["feat_out"])

        regressor_loss, _, _, _, _, _, pred_x_regressor = compute_loss(
            binary_labels,
            orig_labels,
            regressor_output,
            self.regressor_criterion,
            mode="regression",
            pred_binary=pred_binary,
        )
        
        mean_per_col = pred_x_classifier.mean(dim=0)  # shape [1, label_dim]
        # pred_x: shape [N, S]
        
        quantile_25 = pred_x_classifier.quantile(0.25, dim=0)             # shape [S]

        # Masks
        high_mean_mask = (mean_per_col >= 0.2)                 # Rule 1
        mid_mean_mask = (mean_per_col < 0.2) & (mean_per_col > 0.01)  # Rule 2
        low_mean_mask = mean_per_col <= 0.01                   # Rule 3

        # Expand masks to match pred_x shape [N, S]
        high_mask_exp = high_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)
        mid_mask_exp = mid_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)
        low_mask_exp = low_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)

        # Thresholds
        bin_high = (pred_x_classifier > mean_per_col).float()                      # Rule 1 threshold
        bin_mid = (pred_x_classifier > mean_per_col).float()               # Rule 2 threshold
        bin_low = (pred_x_classifier > 0.05).float()                 # Rule 3 all zeros

        # Apply conditional logic
        pred_binary_rounded = torch.where(
            high_mask_exp, bin_high,
            torch.where(
                mid_mask_exp, bin_mid,
                bin_low
            )
        )
        
        # Calculate accuracy using torchmetrics
        #pred_binary_rounded = (pred_x_classifier > mean_per_col).float()
        self.test_accuracy(pred_binary_rounded, binary_labels)

        # Calculate RMSE using torchmetrics
        # If pred_binary < 0.01, set pred_x from regressor as 0, otherwise keep original
        adjusted_pred_x = torch.where(
            pred_binary < 0.01, torch.zeros_like(pred_x_regressor), pred_x_regressor
        )
        # Inverse transform predictions and labels
        pred_x_real = adjusted_pred_x


        # Round pred_x to integers
        rounded_pred_x = torch.round(pred_x_real)
        # Move both to CPU before passing to metric
        rounded_pred_x = rounded_pred_x.cpu()
        orig_labels = orig_labels.cpu()
        self.test_mse(rounded_pred_x, orig_labels)
        
        # Log per-batch loss metrics
        self.log("test/classifier_loss", classifier_loss)
        self.log("test/regressor_loss", regressor_loss)

        # Store predictions and ground truth for CSV export
        self.test_binary_gt.append(binary_labels.cpu().numpy())
        self.test_count_gt.append(orig_labels.cpu().numpy())
        self.test_binary_pred.append(pred_binary_rounded.cpu().numpy())
        self.test_count_pred.append(rounded_pred_x.cpu().numpy())
        

    def on_test_epoch_end(self):
        """
        Save predictions and ground truth to NPZ files at end of test
        """
        # Log accumulated metrics at the end of the epoch
        self.log("test/accuracy", self.test_accuracy.compute())
        self.log("test/rmse", torch.sqrt(self.test_mse.compute()))

        # Reset metrics
        self.test_accuracy.reset()
        self.test_mse.reset()

        # Concatenate collected data
        binary_gt = np.concatenate(self.test_binary_gt)
        count_gt = np.concatenate(self.test_count_gt)
        binary_pred = np.concatenate(self.test_binary_pred)
        count_pred = np.concatenate(self.test_count_pred)

        # Get experiment directory from logger
        if self.logger and hasattr(self.logger, "log_dir"):
            log_dir = self.logger.log_dir
        else:
            log_dir = os.getcwd()

        # Save binary data to NPZ file (first half columns are GT, second half are Pred)
        binary_data = np.hstack([binary_gt, binary_pred])
        binary_npz_path = os.path.join(log_dir, "binary_results.npz")
        np.savez(binary_npz_path, data=binary_data)

        # Save count data to NPZ file (first half columns are GT, second half are Pred)
        count_data = np.hstack([count_gt, count_pred])
        count_npz_path = os.path.join(log_dir, "count_results.npz")
        np.savez(count_npz_path, data=count_data)

        # Clear saved predictions
        self.test_binary_gt = []
        self.test_count_gt = []
        self.test_binary_pred = []
        self.test_count_pred = []

        print(f"Binary results saved to {binary_npz_path}")
        print(f"Count results saved to {count_npz_path}")
    
    def predict_step(self, batch, batch_idx):
        """
        prediction step to acquire predictions of various scenarios
        """
        features, orig_labels, count_labels, binary_labels = batch

        classifier_output = self.classifier(features, orig_labels, binary_labels)
        regressor_output = self.regressor(features, orig_labels, binary_labels)
        
        # From regressor model
        #emb_regressor = self.regressor.fd_x2[0].weight.detach().cpu().numpy().T

        # From classifier model
        #emb_classifier = self.classifier.fd_x2[0].weight.detach().cpu().numpy().T

        # From regressor
        label_emb_regressor = self.regressor.fe0.weight.detach().cpu().numpy().T  # shape: [emb_size, label_dim]

        # From classifier
        label_emb_classifier = self.classifier.fe0.weight.detach().cpu().numpy().T
        
        # From regressor
        feature_regressor_emb = self.regressor.feat_mp_mu.weight.detach().cpu().numpy().T  # shape: [emb_size, label_dim]

        # From classifier
        feature_classifier_emb = self.classifier.feat_mp_mu.weight.detach().cpu().numpy().T



        #label_emb_regressor = self.regressor.fe0.weight.T @ self.regressor.fe0.weight
        #label_emb_classifier = self.classifier.fe0.weight.T @ self.classifier.fe0.weight



        classifier_loss, _, _, _, _, _, pred_x_classifier = compute_loss(
            binary_labels,
            orig_labels,
            classifier_output,
            self.classifier_criterion,
            mode="classification",
        )

        pred_binary = torch.sigmoid(classifier_output["feat_out"])

        regressor_loss, _, _, _, _, _, pred_x_regressor = compute_loss(
                    binary_labels,
                    orig_labels,
                    regressor_output,
                    self.regressor_criterion,
                    mode="regression",
                    pred_binary=pred_binary,
                )
        

        # Calculate accuracy using torchmetrics
        mean_per_col = pred_x_classifier.mean(dim=0)
        #quantile_25 = pred_x_classifier.quantile(0.25, dim=0)             # shape [S]

        # Masks
        high_mean_mask = (mean_per_col >= 0.2)                 # Rule 1
        mid_mean_mask = (mean_per_col < 0.2) & (mean_per_col > 0.01)  # Rule 2
        low_mean_mask = mean_per_col <= 0.01                   # Rule 3

        # Expand masks to match pred_x shape [N, S]
        high_mask_exp = high_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)
        mid_mask_exp = mid_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)
        low_mask_exp = low_mean_mask.unsqueeze(0).expand_as(pred_x_classifier)

        # Thresholds
        bin_high = (pred_x_classifier > mean_per_col).float()                      # Rule 1 threshold
        bin_mid = (pred_x_classifier > mean_per_col).float()               # Rule 2 threshold
        bin_low = (pred_x_classifier > 0.05).float()                     # Rule 3 all zeros

        # Apply conditional logic
        pred_binary_rounded = torch.where(
            high_mask_exp, bin_high,
            torch.where(
                mid_mask_exp, bin_mid,
                bin_low
            )
        )
       
        
        
        #pred_binary_rounded = (pred_x_classifier > mean_per_col).float()
        self.test_accuracy(pred_binary_rounded, binary_labels)

        # Calculate RMSE using torchmetrics
        # If pred_binary < 0.01, set pred_x from regressor as 0, otherwise keep original
        adjusted_pred_x = torch.where(
            pred_binary < 0.01, torch.zeros_like(pred_x_regressor), pred_x_regressor
        )
        # Round pred_x to integers
        # Inverse transform predictions and labels
        pred_x_real = adjusted_pred_x


        # Round pred_x to integers
        rounded_pred_x = torch.round(pred_x_real)
        self.test_mse(rounded_pred_x, orig_labels)

        
        # Store predictions for CSV export
        self.test_binary_gt.append(binary_labels.cpu().numpy())
        self.test_count_gt.append(orig_labels.cpu().numpy())
        self.test_binary_pred.append(pred_x_classifier.cpu().numpy())
        self.test_count_pred.append(rounded_pred_x.cpu().numpy())
        self.test_label_regressor.append(label_emb_regressor)
        self.test_label_classifier.append(label_emb_classifier)
        self.test_feature_classifier.append(feature_classifier_emb)
        self.test_feature_regressor.append(feature_regressor_emb)
 #       # Store locations for later use
        #self.test_loc.append(loc.cpu().numpy())  

    def on_predict_epoch_end(self):
        """
        Save predictions and ground truth to NPZ files at end of test
        """
        
        # Concatenate collected data
        #binary_gt = np.concatenate(self.test_binary_gt)
        #count_gt = np.concatenate(self.test_count_gt)
        binary_pred = np.concatenate(self.test_binary_pred)
        count_pred = np.concatenate(self.test_count_pred)
        label_regressor = np.concatenate(self.test_label_regressor)
        label_classifier = np.concatenate(self.test_label_classifier)
        feature_classifier = np.concatenate(self.test_feature_classifier)
        feature_regressor = np.concatenate(self.test_feature_regressor)
            


        # Get experiment directory from logger
        if self.logger and hasattr(self.logger, "log_dir"):
            log_dir = self.logger.log_dir
        else:
            log_dir = os.getcwd()

        label_regressor_csv_path = os.path.join(log_dir, "label_label_emb_regressor.csv")
        # Assuming label_regressor is a NumPy array
        label_regressor_df = pd.DataFrame(label_regressor)
        label_regressor_df.to_csv(label_regressor_csv_path, index=False)
        label_classifier_csv_path = os.path.join(log_dir, "label_label_emb_classifier.csv")
        label_classifier_df = pd.DataFrame(label_classifier)
        label_classifier_df.to_csv(label_classifier_csv_path, index=False)
        feature_classifier_csv_path = os.path.join(log_dir, "feature_label_emb_classifier.csv")
        feature_classifier_df = pd.DataFrame(feature_classifier)
        feature_classifier_df.to_csv(feature_classifier_csv_path, index=False)
        feature_regressor_csv_path = os.path.join(log_dir, "feature_label_emb_regressor.csv")
        feature_regressor_df = pd.DataFrame(feature_regressor)
        feature_regressor_df.to_csv(feature_regressor_csv_path, index=False)
        # Step 1: Load label_names from the original NPZ
        data_file = np.load('./DATA/ebird_entire/ebird_data0718.npz', allow_pickle=True)
        label_names = data_file["label_names"]  # shape: (num_labels,)
        loc_pred = np.load('./DATA/ebird_entire/loc0718.npz', allow_pickle=True)["loc"]

        # Step 2: Build column names for GT and Pred
        label_names = label_names.tolist()  # convert to regular list if it's not already
        
        # Final column names with longitude and latitude
        column_names = ["longitude", "latitude"] + label_names

        # Save binary data to NPZ file (first half columns are GT, second half are Pred)
        binary_prediction = np.hstack([loc_pred, binary_pred])
        binary_df = pd.DataFrame(binary_prediction, columns=column_names)
        binary_csv_path = os.path.join(log_dir, "presence_probability_results.csv")
        binary_df.to_csv(binary_csv_path, index=False)

        # Save count data to NPZ file (first half columns are GT, second half are Pred)
        # Inverse log1p transform the predicted count data
        count_pred_original_scale = np.expm1(count_pred)
     

        # Round to nearest integer
        count_pred_rounded = np.rint(count_pred_original_scale).astype(int)

        count_prediction = np.hstack([loc_pred, count_pred_rounded])
        count_df = pd.DataFrame(count_prediction, columns=column_names)
        count_csv_path = os.path.join(log_dir, "count_pred_results.csv")
        count_df.to_csv(count_csv_path, index=False)

        # Clear saved predictions
        self.test_binary_gt = []
        self.test_count_gt = []
        self.test_binary_pred = []
        self.test_count_pred = []
        self.test_loc = []

        print(f"Binary results saved to {binary_csv_path}")
        print(f"Count results saved to {count_csv_path}")


    def configure_optimizers(self):
        """
        Configure optimizers for both models
        """
        # Classifier optimizer
        classifier_optimizer = Adam(
            self.classifier.parameters(), lr=float(self.classifier_args.lr)
        )

        classifier_scheduler = StepLR(
            classifier_optimizer,
            step_size=self.classifier_args.lr_decay_times,
            gamma=self.classifier_args.lr_decay,
        )

        # Regressor optimizer
        regressor_optimizer = Adam(
            self.regressor.parameters(), lr=float(self.regressor_args.lr)
        )

        regressor_scheduler = StepLR(
            regressor_optimizer,
            step_size=self.regressor_args.lr_decay_times,
            gamma=self.regressor_args.lr_decay,
        )

        return [
            {
                "optimizer": classifier_optimizer,
                "lr_scheduler": classifier_scheduler,
            },
            {
                "optimizer": regressor_optimizer,
                "lr_scheduler": regressor_scheduler,
            },
        ]
