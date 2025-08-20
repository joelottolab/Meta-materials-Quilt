import os, argparse, yaml, json
def maybe_scale_fast(cfg):
    if cfg.get("FAST_MODE", False): cfg["simulation"]["points"]=max(250,cfg["simulation"]["points"]//10)
    return cfg
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--out",required=True); args=ap.parse_args()
    os.makedirs(args.out, exist_ok=True); cfg=yaml.safe_load(open(args.config)); cfg=maybe_scale_fast(cfg)
    out_cfg=os.path.join(args.out,"resolved_config.json"); json.dump(cfg, open(out_cfg,"w"), indent=2)
    print("Resolved:", out_cfg); print("Next: python3 window_solver.py --config", out_cfg)
if __name__=="__main__": main()
