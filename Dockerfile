FROM ubuntu:24.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ca-certificates \
    git \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/panvar
COPY . .

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPANVAR_BUILD_TESTING=OFF \
    && cmake --build build -j \
    && strip /opt/panvar/build/panvar || true

FROM mambaorg/micromamba:2.3.2

ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER $MAMBA_USER

# Companion tools only: vg (snarls), odgi (inspect graph sort/flip), bcftools (work with the called VCFs),
# and R for the plotting scripts (ggplot2/ggrepel) plus the interactive node-coverage viewer
# (shiny/plotly/data.table/DT). minimap2 is NOT installed (panvar links the minimap2 API statically).
RUN micromamba install --yes --name base --channel conda-forge --channel bioconda \
    vg \
    odgi \
    bcftools \
    r-base \
    r-ggplot2 \
    r-ggrepel \
    r-data.table \
    r-shiny \
    r-plotly \
    r-dt \
    && micromamba clean --all --yes

USER root
COPY --from=builder /opt/panvar/build/panvar /usr/local/bin/panvar
RUN chmod +x /usr/local/bin/panvar

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "panvar"]
CMD ["--help"]
