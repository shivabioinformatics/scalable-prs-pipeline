# Dockerfile for the PRS Pipeline
#
# This image has everything the pipeline needs to run:
#   - Python 3.11
#   - numpy (for the math in calculate_prs.py)
#   - matplotlib (for the visualization plots)
#
# The base Python slim image is tiny (~150MB) and I just pip install
# the two dependencies on top. No heavy bioinformatics tools here
# because the demo uses pure Python scripts. In production you'd
# add bcftools, beagle, crossmap, etc.
#
# Build:  docker build -t prs-pipeline .
# Test:   docker run prs-pipeline python3 -c "import matplotlib; print('OK')"

FROM python:3.11-slim

# I don't want pip caching bloating my image
ENV PIP_NO_CACHE_DIR=1

# install the Python dependencies my scripts need
RUN pip install numpy matplotlib

# copy the pipeline scripts into the image
# nextflow mounts the work dir at runtime, but having the scripts
# baked in means they're always available
COPY bin/ /opt/prs-pipeline/bin/

# make sure the scripts are executable
RUN chmod +x /opt/prs-pipeline/bin/*.py

# default working directory
WORKDIR /data

# label it so I know what this image is for
LABEL maintainer="Shiva Sadeghpour"
LABEL description="PRS Pipeline - Polygenic Risk Score calculation"
LABEL version="1.0"
