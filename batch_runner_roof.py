# MIT License
import os, argparse, yaml, numpy as np, json

def maybe_scale_fast(cfg):
    if cfg.get("FAST_MODE", False):
        cfg["simulation"]["points"] = max(200, cfg["simulation"]["points"]//10)
        cfg["simulation"]["resolution_cells_per_lambda"] = max(10, cfg["simulation"]["resolution_cells_per_lambda"]//10)
        cfg["simulation"]["timesteps"] = max(30000, cfg["simulation"]["timesteps"]//10)
        cfg["simulation"]["pml_thickness_cells"] = max(10, cfg["simulation"]["pml_thickness_cells"]//2)
    return cfg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    cfg = maybe_scale_fast(cfg)
    with open(os.path.join(args.out, "resolved_config.json"), "w") as f:
        json.dump(cfg, f, indent=2)
    print("Resolved config written to", os.path.join(args.out, "resolved_config.json"))
    print("Next: run `python3 stack_reflectance_meep.py --config {}`".format(os.path.join(args.out, "resolved_config.json")))

if __name__ == "__main__":
    main()
