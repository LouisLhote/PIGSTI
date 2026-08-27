"""
Build config/config_hops_custom.txt from the template + pipeline config.

Always sets pathToList from the generated pathogen list.
Applies, in order:
  - hops_parameters: dict of key=value overrides (matches config_hops2.0.txt)
  - hops_malt_index  -> index=
  - hops_resources   -> resources=
(pathogen list wins for pathToList; top-level hops_malt_index / hops_resources
override the same keys if also present under hops_parameters.)
"""
import os
import sys
import yaml

original = sys.argv[1]
new_config = sys.argv[2]
new_pathogen_list = sys.argv[3]
config_yaml = sys.argv[4] if len(sys.argv) > 4 else "config/config.yaml"

cfg = {}
if config_yaml and os.path.exists(config_yaml):
    try:
        with open(config_yaml, "r", encoding="utf-8") as cf:
            cfg = yaml.safe_load(cf) or {}
    except Exception:
        cfg = {}

overrides = {}
hp = cfg.get("hops_parameters")
if isinstance(hp, dict):
    for k, v in hp.items():
        if v is None or v == "":
            continue
        overrides[str(k)] = str(v).strip()

if cfg.get("hops_malt_index"):
    overrides["index"] = str(cfg["hops_malt_index"]).strip()
if cfg.get("hops_resources"):
    overrides["resources"] = str(cfg["hops_resources"]).strip()

if cfg.get("hops_parallel"):
    per_job = cfg.get("hops_threads_per_job")
    if per_job is not None and str(per_job).strip() != "":
        overrides["threadsMalt"] = str(int(per_job))
    else:
        jobs = max(1, int(cfg.get("hops_parallel_jobs", 2)))
        base_threads = int(overrides.get("threadsMalt", 15))
        overrides["threadsMalt"] = str(max(1, base_threads // jobs))
    per_job_mem = cfg.get("hops_max_memory_malt_per_job")
    if per_job_mem is not None and str(per_job_mem).strip() != "":
        overrides["maxMemoryMalt"] = str(int(per_job_mem))
    elif "maxMemoryMalt" in overrides:
        jobs = max(1, int(cfg.get("hops_parallel_jobs", 2)))
        base_mem = int(overrides["maxMemoryMalt"])
        overrides["maxMemoryMalt"] = str(max(64, base_mem // jobs))

if cfg.get("hops_malt_mmap"):
    overrides["memoryMode"] = "map"

# HOPS requires Java -Xmx (hops_heap_gb) >= maxMemoryMalt; leave JVM headroom.
heap_gb = max(1, int(cfg.get("hops_heap_gb", 800)))
malt_heap_cap = max(64, heap_gb - 20)
if "maxMemoryMalt" in overrides:
    malt_mem = int(overrides["maxMemoryMalt"])
    if malt_mem > malt_heap_cap:
        print(
            f"[create_hops_config] capping maxMemoryMalt {malt_mem} -> {malt_heap_cap} "
            f"(hops_heap_gb={heap_gb}; HOPS needs -Xmx >= maxMemoryMalt)",
            file=sys.stderr,
        )
        overrides["maxMemoryMalt"] = str(malt_heap_cap)

overrides["pathToList"] = new_pathogen_list

with open(original, encoding="utf-8") as f:
    lines = f.readlines()

out_lines = []
for line in lines:
    raw = line.rstrip("\n")
    stripped = raw.lstrip()
    if not stripped or stripped.startswith("#"):
        out_lines.append(line)
        continue
    if "=" not in raw:
        out_lines.append(line)
        continue
    key = raw.split("=", 1)[0].strip()
    if key in overrides:
        out_lines.append(f"{key}={overrides[key]}\n")
    else:
        out_lines.append(line)

os.makedirs(os.path.dirname(os.path.abspath(new_config)) or ".", exist_ok=True)
with open(new_config, "w", encoding="utf-8") as f:
    f.writelines(out_lines)
