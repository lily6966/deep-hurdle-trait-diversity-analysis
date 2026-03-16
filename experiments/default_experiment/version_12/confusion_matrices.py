
exp.dir = "lightning_logs/low/version_0"
binary_path = os.path.join(exp_dir, "binary_results.npz")
binary_data = np.load(binary_path)["data"]