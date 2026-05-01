FROM hkubal/clair3:v2.0.1_gpu

RUN mamba install -n clair3_v2 -c bioconda -c conda-forge -y \
        bedtools=2.31.1 \
        ucsc-bedgraphtobigwig=482 \
    && mamba clean --all -y
