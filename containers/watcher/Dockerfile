FROM hkubal/clair3:v1.2.0_gpu

RUN mamba install -n clair3 -c bioconda -c conda-forge -y \
        ucsc-bigwigmerge=482 \
        ucsc-bedgraphtobigwig=482 \
    && mamba clean --all -y
