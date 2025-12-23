FROM python:3.11-slim

WORKDIR /opt/ensembl

COPY ../../../pyproject.toml .
RUN pip install --no-cache-dir .

COPY src/python /opt/ensembl/src/python
ENV PYTHONPATH=/opt/ensembl/src/python
