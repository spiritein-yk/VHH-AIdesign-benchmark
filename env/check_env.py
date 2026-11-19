#!/usr/bin/env python3
"""Check configured executables, inputs and imports without running a model."""
import argparse
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime import STAGE_PATHS, location

FILES = {
    "ppiflow": [("repo", "tools/PPIFlow/sample_antibody_nanobody.py")],
    "iggm": [("repo", "design.py")],
    "proteinmpnn": [("repo", "protein_mpnn_run.py"), ("weights", "v_48_020.pt")],
    "solublempnn": [("repo", "protein_mpnn_run.py"), ("weights", "v_48_020.pt")],
    "antibmpnn": [("repo", "Running_AntiBMPNN_run.py")],
    "packing": [("repo", "sc_utils.py"), ("repo", "model_params/ligandmpnn_sc_v_32_002_16.pt")],
    "silent": [("tools", name) for name in ("silentfrompdbs", "silentslice", "silentls", "silentrename")],
    "af2": [("repo", "alphafold/model/model.py"), ("silent_tools", "silent_tools.py"),
            ("weights", "params/params_model_3_multimer_v2.npz")],
    "af3": [("repo", "run_alphafold.py")],
}
PROBES = {
    "ppiflow": "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())",
    "iggm": "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())",
    "pdbfixer": "import openmm, pdbfixer; print(openmm.__version__)",
    "proteinmpnn": "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())",
    "solublempnn": "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())",
    "antibmpnn": "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())",
    "packing": ("import numpy, torch; print(torch.__version__, torch.version.cuda); "
                "assert torch.cuda.is_available(), 'CUDA GPU is required for packing'"),
    "silent": "import numpy; print(numpy.__version__)",
    "af2": ("import alphafold, jax, pyrosetta, silent_tools; print(jax.__version__, jax.devices()); "
            "assert any(d.platform == 'gpu' for d in jax.devices()), 'JAX GPU device unavailable'"),
    "af3": ("import alphafold3, jax; print(jax.__version__, jax.devices()); "
            "assert any(d.platform == 'gpu' for d in jax.devices()), 'JAX GPU device unavailable'"),
}


def check(stage, files_only):
    paths = {key: location(stage, key) for key in STAGE_PATHS[stage]}
    for key, filename in FILES.get(stage, ()):
        file = paths[key] / filename
        if not file.is_file() or not os.access(file, os.R_OK):
            raise ValueError(f"[{stage}].{key} is missing readable {filename}")
    if not files_only and stage in PROBES:
        env = dict(os.environ, PYTHONNOUSERSITE="1")
        if stage == "af2":
            env["PYTHONPATH"] = os.pathsep.join(str(paths[key]) for key in ("repo", "silent_tools"))
        result = subprocess.run([paths["python"], "-c", PROBES[stage]], cwd=paths.get("repo"),
                                env=env, check=True, timeout=60, text=True, capture_output=True)
        print(f"{stage} runtime  {result.stdout.strip()}")
    print(f"{stage} PASS")


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", action="append", choices=STAGE_PATHS, required=True)
    parser.add_argument("--files-only", action="store_true")
    args = parser.parse_args()
    failed = False
    for stage in args.stage:
        try:
            check(stage, args.files_only)
        except (ValueError, OSError, subprocess.SubprocessError) as exc:
            details = getattr(exc, "stderr", None) or str(exc)
            print(f"{stage} FAIL  {details.strip()}", file=sys.stderr)
            failed = True
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    run()
