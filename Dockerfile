# Pinned Bioconductor release: R + Bioconductor versions are fixed, so the
# analysis gives the same result on any machine (laptop, HPC, CI).
FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL org.opencontainers.image.source="https://github.com/NicoCosta4/methylation-arrays" \
      org.opencontainers.image.description="Reproducible Illumina EPIC methylation array analysis" \
      org.opencontainers.image.licenses="MIT"

COPY docker/install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm /tmp/install_packages.R

WORKDIR /project
CMD ["bash"]
