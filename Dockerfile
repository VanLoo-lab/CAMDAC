# Use a lightweight R base image
FROM rocker/r-ver:4.4.0
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
    && apt-get install -y default-jre git curl bzip2 pandoc

# Install alleleCounter
RUN git clone https://github.com/cancerit/alleleCount.git \
    && cd alleleCount && bash ./setup.sh /usr/local
RUN echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
RUN echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bashrc

# Set WORKDIR
WORKDIR /app

# Install CAMDAC
RUN R -q -e 'install.packages("remotes")'
RUN R -q -e 'install.packages("devtools")'

# Create ~/.R directory and configure Makevars file for R-specific flags.
#   Required for bioconductor packages used in Battenberg install.
RUN mkdir -p /root/.R && \
    echo 'CXXFLAGS=-Wall -Wno-format-security' >> /root/.R/Makevars && \
    echo 'CFLAGS=-Wall -Wno-format-security' >> /root/.R/Makevars

# Copy only DESCRIPTION and (optionally) NAMESPACE to install deps first
COPY DESCRIPTION NAMESPACE* ./
RUN R -q -e 'remotes::install_deps(".", dependencies = TRUE, upgrade = "never", lib = "/usr/local/lib/R/site-library")'

# Now copy the rest of your project (code layer)
# COPY . /app

# Install CAMDAC References from repository [Deprecated]
# Optionally install CAMDAC if it's not a local package (e.g., separate install)
# RUN Rscript -e 'remotes::install_github("VanLoo-lab/CAMDAC@wgbs", lib="/usr/local/lib/R/site-library")'
#RUN R -q -e 'library(CAMDAC);CAMDAC::download_pipeline_files(bsseq="wgbs", directory="pipeline_files/")'

# Set command to be use
CMD ["/usr/local/bin/Rscript"]
