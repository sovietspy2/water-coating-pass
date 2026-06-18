FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y git bash ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN git clone https://github.com/sovietspy2/water-coating-pass.git

WORKDIR /opt/water-coating-pass

RUN chmod +x installer/install.sh testing/test.sh && \
    cd installer && \
    ./install.sh

CMD ["/bin/bash", "-lc", "mkdir -p /data/test_runs && cd /opt/water-coating-pass/testing && ./test.sh -p /data/test_runs"]