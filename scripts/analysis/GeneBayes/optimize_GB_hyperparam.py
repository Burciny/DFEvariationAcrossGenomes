import optuna
import subprocess
import re
import argparse
import optuna.visualization as vis

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--likelihood", required=True)
    parser.add_argument("--feature", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--n_trials", type=int, default=100, help="Number of Optuna trials")

    return parser.parse_args()


args = parse_args()
SCRIPT = "genebayes_shet_modified_filter.py" # Change it for the species if there are script differences

def extract_val_loss(output):
    """
    Extract validation loss from standard output.
      1. Best-iteration summary line
      2. Otherwise: take min of all val_loss occurrences (should i???)
    """

    # 1) Look for final best-iteration summary
    best_iter_match = re.search(r"Best iteration.*?val_loss=([0-9.]+)", output)
    if best_iter_match:
        return float(best_iter_match.group(1))

    # 2) Collect all val_loss values printed during training
    all_losses = re.findall(r"val_loss=([0-9.]+)", output)

    if not all_losses:
        print("No val_loss found in output. Full output:")
        print(output)
        raise RuntimeError("Could not extract validation loss.")

    # Return minimum — best validation performance
    return min(float(v) for v in all_losses)


def run_model(params):

    cmd = [
        "python", SCRIPT,
        "--response", args.likelihood,
        "--features", args.feature,
        "--out", f"{args.prefix}_optuna_temp",
        "--n_jobs", "8",
        "--batch_size", "100",
        "--n_integration_pts", "10000",
        "--lr", str(params["lr"]),
        "--max_depth", str(params["max_depth"]),
        "--min_child_weight", str(params["min_child_weight"]),
        "--subsample", str(params["subsample"]),
        "--n_trees_per_iteration", str(params["n_trees_per_iteration"]),
        "--reg_alpha", str(params["reg_alpha"]),
        "--reg_lambda", str(params["reg_lambda"])
    ]

    # capture output
    output = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)

    # this handles both best-iteration and fallback-to-min cases
    loss = extract_val_loss(output)

    return loss

### Narrow the ranges depending on the previous best grid search run,
def objective(trial):
    """Optuna objective function for Bayesian optimization. Suggestions change depending on species parameters"""
    params = {
        "lr": trial.suggest_float("lr", 0.01, 0.1, log=True), 
        "max_depth": trial.suggest_int("max_depth", 2, 4),
        "min_child_weight": trial.suggest_float("min_child_weight", 0.8, 5, log=True), 
        "subsample": trial.suggest_float("subsample", 0.6, 0.9), 
        "n_trees_per_iteration": trial.suggest_int("n_trees_per_iteration", 1, 4), 
        "reg_alpha": trial.suggest_float("reg_alpha",  1e-8, 5, log=True), 
        "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 5, log=True), 
    }

    return run_model(params)


if __name__ == "__main__":
    study = optuna.create_study(direction="minimize")
"""Start parameters change depending on grid search"""
    enqueue_params = {
        "lr": 0.01,
        "max_depth": 3,
        "min_child_weight": 1,
        "subsample": 0.8,
        "n_trees_per_iteration": 1,
        "reg_alpha": 1,
        "reg_lambda": 1,
    }
    study.enqueue_trial(enqueue_params)
    study.optimize(objective, n_trials=args.n_trials)

    print("\n=== Best hyperparameters ===")
    print(study.best_params)
    print("Best validation loss:", study.best_value)

    # final training with best hyperparameters
    best = study.best_params

    print("\nTraining final model on best params...")
    final_cmd = [
        "python", SCRIPT,
        "--response", args.likelihood,
        "--features", args.feature,
        "--out", f"{args.prefix}_best_model",
        "--n_jobs", "8",
        "--batch_size", "100",
        "--n_integration_pts", "10000",
        "--early_stopping_iter", "10",
        "--lr", str(best["lr"]),
        "--max_depth", str(best["max_depth"]),
        "--min_child_weight", str(best["min_child_weight"]),
        "--subsample", str(best["subsample"]),
        "--n_trees_per_iteration", str(best["n_trees_per_iteration"]),
        "--reg_alpha", str(best["reg_alpha"]),
        "--reg_lambda", str(best["reg_lambda"])
    ]

    subprocess.check_call(final_cmd)
    print("Done.")

### Optimization summary plots ####
#    fig1 = vis.plot_optimization_history(study)
#    fig1.write_image(f"{args.prefix}_opt_history.png")

    # 2. Hyperparameter importance
#    fig2 = vis.plot_param_importances(study)
#    fig2.write_image(f"{args.prefix}_hyperparam_importance.png")

    # 3. Parallel coordinate plot
#    fig3 = vis.plot_parallel_coordinate(study)
#    fig3.write_image(f"{args.prefix}_parallel_coord.png")


