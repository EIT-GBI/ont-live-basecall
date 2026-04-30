FROM hkubal/clair3:v2.0.1_gpu

RUN mamba install -n clair3_v2 -c bioconda -c conda-forge -y \
        ucsc-bigwigmerge=482 \
        ucsc-bedgraphtobigwig=482 \
    && mamba clean --all -y
