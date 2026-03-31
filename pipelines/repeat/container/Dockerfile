FROM python:3.11-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Clone and build
WORKDIR /build
RUN git clone --depth 1 --branch dev/python_repo https://github.com/Ensembl/ensembl-anno.git

# Build the wheel - cd into the cloned repository
WORKDIR /build/ensembl-anno

# Upgrade pip and install build tools
RUN pip install --no-cache-dir --upgrade pip setuptools wheel
RUN pip install --no-cache-dir .

LABEL maintainer="Ensembl <ensembl@ebi.ac.uk.com>"
LABEL description="Optimized Docker image with ensembl-anno for RepeatMasker"

# Set the working directory back to root
WORKDIR /

# Verify installation
RUN python -c "from importlib.metadata import distribution; pkg = distribution('ensembl-anno'); print(f'✓ Package installed: {pkg.metadata[\"Name\"]} v{pkg.version}')"

# Verify the correct import path
RUN python -c "import ensembl.tools.anno; print('✓ Successfully imported ensembl.tools.anno')"

# Show available submodules
RUN python -c "import ensembl.tools.anno; import os; anno_path = ensembl.tools.anno.__path__[0]; print(f'Available annotation modules:'); [print(f'  - {d}') for d in os.listdir(anno_path) if os.path.isdir(os.path.join(anno_path, d)) and not d.startswith('_')]"

# Default command
CMD ["/bin/bash"]