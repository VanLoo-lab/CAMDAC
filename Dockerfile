# Use a lightweight R base image
FROM rocker/r-ver:4.3.0
# Set the working directory
WORKDIR /app

# Install system dependencies (if needed)
# For example, if your package needs external libraries, install them here
RUN chmod 1777 /tmp
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    && apt-get install -y libxml2 libodbc1 \
    && apt-get install -y libz-dev libbz2-dev bzip2-doc zlib1g-dev \
    && apt-get install -y liblzma-dev libcurl4-openssl-dev wget \
    && apt-get install -y default-jre git curl bzip2

# Install alleleCounter
RUN git clone https://github.com/cancerit/alleleCount.git \
    && cd alleleCount && bash ./setup.sh /usr/local
RUN echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
RUN echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bashrc

# Install CAMDAC
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'remotes::install_github("VanLoo-lab/CAMDAC@wgbs", lib="/usr/local/lib/R/site-library")'

# Install CAMDAC References from repository [Deprecated]
#RUN R -q -e 'library(CAMDAC);CAMDAC::download_pipeline_files(bsseq="wgbs", directory="pipeline_files/")'

# Install CAMDAC references from local WGBS files
# Copy local WGBS pipeline files and extract for CAMDAC
RUN mkdir -p /opt/pipeline_files
COPY ./camdac_wgbs_pipeline_files.tar.gz /opt/camdac_wgbs_pipeline_files.tar.gz
RUN tar -zxvf /opt/camdac_wgbs_pipeline_files.tar.gz -C /opt/pipeline_files
RUN rm /opt/camdac_wgbs_pipeline_files.tar.gz

# Set the working directory
WORKDIR /app

# Set command to be use
CMD ["/usr/local/bin/Rscript"]
