"""Paths and process invocation for the release entry points."""
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TARGETS = ("PDGFRB", "IL7RA")
FLOWS = ("boltzgen", "ppiflow", "rfantibody", "iggm")
MODELS = ("proteinmpnn", "solublempnn", "antibmpnn")

# Each dependency specifies its environment variable, path type and fallback.
STAGE_PATHS = {
    "boltzgen": {"exe": ("BOLTZGEN_EXE", "exe", "boltzgen"),
                 "cache": ("BOLTZGEN_CACHE", "dir", None)},
    "ppiflow": {"python": ("PPIFLOW_PYTHON", "exe", sys.executable),
                "repo": ("PPIFLOW_REPO", "dir", None),
                "weights": ("PPIFLOW_WEIGHTS", "file", None)},
    "rfantibody": {"exe": ("RFANTIBODY_EXE", "exe", "rfdiffusion"),
                   "repo": ("RFANTIBODY_REPO", "dir", None)},
    "iggm": {"python": ("IGGM_PYTHON", "exe", sys.executable),
             "repo": ("IGGM_REPO", "dir", None)},
    "pdbfixer": {"python": ("PDBFIXER_PYTHON", "exe", sys.executable)},
    "proteinmpnn": {"python": ("PROTEINMPNN_PYTHON", "exe", sys.executable),
                    "repo": ("PROTEINMPNN_REPO", "dir", None),
                    "weights": ("PROTEINMPNN_WEIGHTS", "dir", None)},
    "solublempnn": {"python": ("PROTEINMPNN_PYTHON", "exe", sys.executable),
                    "repo": ("PROTEINMPNN_REPO", "dir", None),
                    "weights": ("SOLUBLEMPNN_WEIGHTS", "dir", None)},
    "antibmpnn": {"python": ("ANTIBMPNN_PYTHON", "exe", sys.executable),
                  "repo": ("ANTIBMPNN_REPO", "dir", None)},
    "packing": {"python": ("LIGANDMPNN_PYTHON", "exe", sys.executable),
                "repo": ("LIGANDMPNN_REPO", "dir", None)},
    "silent": {"python": ("SILENT_PYTHON", "exe", sys.executable),
               "tools": ("SILENT_TOOLS_DIR", "dir", None)},
    "rosetta": {"exe": ("ROSETTA_EXE", "exe", "rosetta_scripts.mpi.linuxgccrelease"),
                "mpirun": ("MPI_RUN", "exe", "mpirun"),
                "dalphaball": ("DALPHABALL_EXE", "exe", None),
                "psipred": ("PSIPRED_EXE", "exe", None)},
    "af2": {"python": ("AF2_PYTHON", "exe", sys.executable),
            "repo": ("AF2_REPO", "dir", None),
            "silent_tools": ("SILENT_TOOLS_DIR", "dir", None),
            "weights": ("AF2_DATA_DIR", "dir", None)},
    "af3": {"python": ("AF3_PYTHON", "exe", sys.executable),
            "repo": ("AF3_REPO", "dir", None),
            "weights": ("AF3_WEIGHTS", "dir", None),
            "databases": ("AF3_DATABASES", "dir", None),
            "pdb_database": ("AF3_PDB_DATABASE", "file", None)},
}


def location(stage, key, *, check=True):
    variable, kind, default = STAGE_PATHS[stage][key]
    raw = os.environ.get(variable) or default
    if raw is None and key == "weights" and stage in ("proteinmpnn", "solublempnn"):
        subdirectory = "soluble_model_weights" if stage == "solublempnn" else "vanilla_model_weights"
        raw = str(location(stage, "repo", check=check) / subdirectory)
    if raw is None:
        if not check:
            return Path(f"${variable}")
        raise ValueError(f"Set {variable} before running {stage}")

    raw = os.path.expanduser(raw)
    if kind == "exe":
        executable = shutil.which(raw)
        if executable:
            return Path(os.path.abspath(executable))
        if not check:
            return Path(raw)
        raise ValueError(f"{variable} is not executable or on PATH: {raw}")

    if not check and raw.startswith("$"):
        return Path(raw)
    path = Path(os.path.abspath(raw))
    valid = path.is_dir() if kind == "dir" else path.is_file()
    if check and (not valid or not os.access(path, os.R_OK)):
        raise ValueError(f"{variable} is not a readable {kind}: {path}")
    return path


def invoke(command, *, cwd=None, dry_run=False, env=None):
    command = [str(part) for part in command]
    print(shlex.join(command), flush=True)
    if not dry_run:
        subprocess.run(command, cwd=cwd, env=env, check=True)


def main(function):
    try:
        function()
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        sys.exit(str(exc))
