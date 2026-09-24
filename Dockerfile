FROM rust:1-bookworm AS rust-builder

WORKDIR /build
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip python3-venv \
    && python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir maturin \
    && rm -rf /var/lib/apt/lists/*
ENV PATH="/opt/venv/bin:${PATH}"

COPY Cargo.toml pyproject.toml ./
COPY src ./src
COPY dna_compare ./dna_compare

RUN maturin build --release -o /wheels

FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY --from=rust-builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/*.whl && rm -rf /wheels

COPY . .

EXPOSE 8765

CMD ["python", "main.py", "serve", "--host", "0.0.0.0", "--port", "8765", "--no-open"]
