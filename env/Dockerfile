# ═══════════════════════════════════════════════════════════════════════════
#  CaptureForge — Docker image (Mac/CI development)
#  Build:  docker build -t captureforge:0.1.0 -f env/Dockerfile .
#  Uses micromamba for a small, reproducible, pinned conda environment.
# ═══════════════════════════════════════════════════════════════════════════
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="CaptureForge" \
      org.opencontainers.image.description="Hybridization-capture probe design pipeline" \
      org.opencontainers.image.version="0.1.0" \
      org.opencontainers.image.licenses="MIT"

USER root
# linux-64 builds are native here (no emulation), so the same yml resolves.
COPY env/captureforge.yml /tmp/captureforge.yml
RUN micromamba install -y -n base -f /tmp/captureforge.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Sanity-check the toolchain is on PATH at build time.
RUN bedtools --version && samtools --version | head -1 && bowtie2 --version | head -1 && \
    genmap --version && dustmasker -version && primer3_core --help 2>&1 | head -1 && \
    python -c "import pysam, Bio, pyfaidx, primer3, pandas, openpyxl, jinja2, matplotlib; print('py libs OK')"

WORKDIR /work
CMD ["/bin/bash"]
