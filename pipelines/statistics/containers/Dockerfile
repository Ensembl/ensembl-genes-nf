# Use official slim Python 3.11 image
FROM python:3.11-slim

# Set working directory
WORKDIR /opt/ensembl

# Copy requirements and install dependencies
COPY requirements.txt .
RUN pip install --upgrade pip setuptools wheel \
    && pip install -r requirements.txt


# Default to bash so scripts can be run externally
CMD ["bash"]