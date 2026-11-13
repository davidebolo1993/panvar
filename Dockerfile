# panvar, with every tool its scripts and dev/ checks call.

# Built on the runtime image, so the binary matches its C and C++ libraries.
FROM mambaorg/micromamba:2.3.2 AS builder
USER root
# bzip2, lzma and libdeflate give the vendored htslib its CRAM codecs and faster BGZF decompression
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake zlib1g-dev \
        libbz2-dev liblzma-dev libdeflate-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/panvar
COPY . .
# Only strip may fail; a configure or compile error stops the build.
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPANVAR_BUILD_TESTING=OFF \
    && cmake --build build -j \
    && { strip build/panvar || true; }

FROM mambaorg/micromamba:2.3.2
USER root
# column (bsdextrautils) for the scripts' tables; the libraries panvar links for CRAM codecs
RUN apt-get update && apt-get install -y --no-install-recommends bsdextrautils libbz2-1.0 liblzma5 libdeflate0 \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
#   graphs and VCFs   vg (snarls for bubble --snarls), odgi, bcftools, htslib (bgzip, tabix)
#   dev/ checks       minimap2 with k8 and paftools.js, samtools, svim-asm, truvari, gemma; numpy
#   recruit dev/      read simulators (wgsim, ART, pbsim3), BBDuk (bbmap), Locityper with jellyfish
#                     (mirabait, from MIRA, is left out: dev/recruit_baselines.sh skips it when absent)
#   R                 ggplot2, ggrepel (plots); data.table, shiny, plotly (viewer); SKAT, SPAtest,
#                     susieR (dev/validate_associate)
RUN micromamba install --yes --name base --channel conda-forge --channel bioconda \
    python numpy \
    vg odgi bcftools htslib \
    minimap2 k8 samtools svim-asm truvari gemma \
    wgsim art pbsim3 bbmap locityper kmer-jellyfish \
    r-base r-ggplot2 r-ggrepel r-data.table r-shiny r-plotly r-skat r-spatest r-susier \
    && micromamba clean --all --yes

USER root
COPY --from=builder /opt/panvar/build/panvar /usr/local/bin/panvar
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "panvar"]
CMD ["--help"]
