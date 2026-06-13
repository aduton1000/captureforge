-- Lmod modulefile stub for CaptureForge (edit INSTALL_ROOT for your cluster)
help([[
CaptureForge 0.1.0 — hybridization-capture probe (bait) design pipeline.
Loads Nextflow + Apptainer on PATH and exposes `captureforge-run`.
Run:  captureforge-run -params-file my_panel.params.yaml --slurm_partition <queue>
]])

whatis("Name: CaptureForge")
whatis("Version: 0.1.0")
whatis("Description: General-purpose hybridization-capture probe design pipeline (Nextflow)")

local install_root = "/shared/apps/captureforge"   -- EDIT for your cluster
local refs_root    = "/shared/refs/captureforge/GRCh38"

-- Dependencies (adjust to your cluster's module names)
depends_on("nextflow")
depends_on("apptainer")

setenv("CAPTUREFORGE_HOME", install_root)
setenv("CAPTUREFORGE_REFS", refs_root)
setenv("CAPTUREFORGE_SIF", pathJoin(install_root, "captureforge.sif"))
setenv("CAPTUREFORGE_APPTAINER_CACHE", pathJoin(os.getenv("HOME") or "/tmp", ".apptainer/captureforge"))

prepend_path("PATH", pathJoin(install_root, "share/bin"))
prepend_path("PATH", pathJoin(install_root, "bin"))
