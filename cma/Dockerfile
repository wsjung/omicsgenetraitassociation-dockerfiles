FROM thomaschln/r-devtools

WORKDIR /app
ADD install_packages.R /app
RUN Rscript /app/install_packages.R

COPY cma.R /app
