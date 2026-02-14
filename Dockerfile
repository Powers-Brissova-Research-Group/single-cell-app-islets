FROM rocker/shiny:4.4.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libhdf5-dev \
    libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('shiny', 'shinyjs', 'Seurat', 'ggplot2', 'Matrix', \
    'dplyr', 'shinydashboard', 'shinyWidgets', 'DT', 'shinycssloaders'), repos='https://cloud.r-project.org/')"

# Copy app files
WORKDIR /srv/shiny-server
COPY app.R .
COPY DATA/ ./DATA/
COPY www/ ./www/

# Expose port
EXPOSE 3838

# Run app on all interfaces
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=3838, launch.browser = FALSE)"]
