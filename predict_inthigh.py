import argparse
import pytorch_lightning as pl
from lightning_models0617 import DeepHurdleModule
from utils import load_config, load_data, get_dataloaders, get_dataloaders_esrd, load_esrd_inthigh
import subprocess
import os


def main(args):
    # If checkpoint is provided, load model directly without config
    if args.ckpt:
        model = DeepHurdleModule.load_from_checkpoint(args.ckpt)

        # Set up basic trainer for testing
        trainer = pl.Trainer(
            accelerator="auto",
            deterministic=True,
        )

        # Load data for testing
        config =  model.hparams["config"]
        data = load_data(config)
        esrd = load_esrd_inthigh()
        test_dataloader = get_dataloaders(data, config, split="test")
        esrd_dataloader = get_dataloaders_esrd(esrd, config)

        # Test model
        trainer.test(model, dataloaders=test_dataloader)
        
         
        # Run analysis on test results
        # Default experiment directory is the parent directory of checkpoint
        exp_dir = trainer.logger.log_dir
        print(f"Experiment directory: {exp_dir}")
        run_analysis(exp_dir)
      
            # Make pridictions on ESRD data
        trainer.predict(model, dataloaders=esrd_dataloader)

    else:
        # Regular training flow
        config = load_config(args.config)
       
        data = load_data(config)
        train_dataloader = get_dataloaders(data, config, split="train")
        val_dataloader = get_dataloaders(data, config, split="val")
        test_dataloader = get_dataloaders(data, config, split="test")
        
        pl.seed_everything(config["data"]["random_seed"])

        # Initialize model
        model = DeepHurdleModule(config)

        # Set up experiment logging
        exp_name = args.exp if args.exp else "default_experiment"
        logger = pl.loggers.TensorBoardLogger("experiments", name=exp_name)
        exp_dir = os.path.join("experiments", exp_name, f"version_{logger.version}")

        trainer = pl.Trainer(
            max_epochs=config["training"]["max_epoch"] * 2,
            accelerator="auto",
            deterministic=True,
            logger=logger,
        )

        # Train model
        trainer.fit(
            model, train_dataloaders=train_dataloader, val_dataloaders=val_dataloader
        )
        # Test model
        trainer.test(model, dataloaders=test_dataloader)


        # Run analysis on results
        run_analysis(exp_dir)


def run_analysis(exp_dir):
    """Run the analysis script on the experiment results."""
    try:
        print(f"Running analysis on results in {exp_dir}...")
        # Capture both stdout and stderr
        result = subprocess.run(
            ["python", "analyze_results.py", "--exp_dir", exp_dir],
            check=True,
            capture_output=True,
            text=True,
        )

        # Save the output to a log file in the experiment directory
        log_path = os.path.join(exp_dir, "analysis_log.txt")
        with open(log_path, "w") as log_file:
            log_file.write(result.stdout)
            if result.stderr:
                log_file.write("\n\nERRORS/WARNINGS:\n")
                log_file.write(result.stderr)

        # Still print to console for visibility
        print(result.stdout)
        if result.stderr:
            print("ERRORS/WARNINGS:")
            print(result.stderr)

        print(f"Analysis completed successfully. Output saved to {log_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error running analysis: {e}")
        # Save error output if available
        log_path = os.path.join(exp_dir, "analysis_error_log.txt")
        with open(log_path, "w") as log_file:
            log_file.write(f"Analysis failed with error code: {e.returncode}\n\n")
            if e.stdout:
                log_file.write("STDOUT:\n")
                log_file.write(e.stdout)
            if e.stderr:
                log_file.write("\n\nSTDERR:\n")
                log_file.write(e.stderr)
        print(f"Error details saved to {log_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Deep Hurdle Population Model")
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to config file (not needed if checkpoint is provided)",
    )
    parser.add_argument("--exp", type=str, help="Name for the experiment")
    parser.add_argument(
        "--ckpt", type=str, help="Path to checkpoint file for testing only"
    )
    args = parser.parse_args()

    main(args)
