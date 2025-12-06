FROM python:3.10-slim

LABEL description="Container for unique reads merger with Python 3.10 and polars" \
      maintainer="EnsEMBL Genebuild"

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    polars==0.20.0

# Set working directory
WORKDIR /work

# Default command
CMD ["/bin/bash"]
