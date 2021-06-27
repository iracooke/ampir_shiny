FROM rocker/shiny-verse

MAINTAINER Ira Cooke "ira.cooke@jcu.edu.au"

RUN R -e "install.packages(c('ampir'), repos='https://cloud.r-project.org/')"

# copy the app to the image
RUN mkdir /srv/shiny-server/ampir

COPY ui.R /srv/shiny-server/ampir
COPY server.R /srv/shiny-server/ampir
COPY about.md /srv/shiny-server/ampir
COPY ampir_hex.png /srv/shiny-server/ampir
COPY google-analytics.html /srv/shiny-server/ampir
COPY instructions.html /srv/shiny-server/ampir
COPY style.css /srv/shiny-server/ampir

COPY shiny.conf /etc/shiny-server/shiny-server.conf


