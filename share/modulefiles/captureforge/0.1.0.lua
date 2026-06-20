-- Lmod modulefile stub for CaptureForge (edit INSTALL_ROOT for your cluster).
-- OPTIONAL: only for sites that have Lmod. On a no-module cluster you do NOT need
-- this — put <install_root>/share/bin on PATH and run `captureforge-run` directly.
help([[
CaptureForge 0.1.0 — hybridization-capture probe (bait) design pipeline.
Loads Nextflow + a container engine (Apptainer/Singularity) on PATH and exposes `captureforge-run`.
Run:  captureforge-run -params-file my_panel.params.yaml --slurm_partition <queue>
]])

whatis("Name: CaptureForge")
whatis("Version: 0.1.0")
whatis("Description: General-purpose hybridization-capture probe design pipeline (Nextflow)")

local install_root = "/shared/apps/captureforge"   -- EDIT for your cluster
local refs_root    = "/shared/refs/captureforge/GRCh38"

-- Dependencies (adjust to your cluster's module names / container engine)
depends_on("nextflow")
depends_on("apptainer")   -- on a SingularityCE site use: depends_on("singularity")

setenv("CAPTUREFORGE_HOME", install_root)
setenv("CAPTUREFORGE_REFS", refs_root)
setenv("CAPTUREFORGE_SIF", pathJoin(install_root, "captureforge.sif"))
-- Engine the captureforge-run wrapper composes into -profile hpc_slurm,<engine>.
-- The wrapper defaults to singularity; set apptainer here on an Apptainer site.
setenv("CAPTUREFORGE_ENGINE", "apptainer")
setenv("CAPTUREFORGE_APPTAINER_CACHE", pathJoin(os.getenv("HOME") or "/tmp", ".apptainer/captureforge"))

prepend_path("PATH", pathJoin(install_root, "share/bin"))
prepend_path("PATH", pathJoin(install_root, "bin"))
