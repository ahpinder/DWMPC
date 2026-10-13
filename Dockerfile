# syntax=docker/dockerfile:1
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates python3 python3-dev python3-pip \
        libeigen3-dev libyaml-cpp-dev pybind11-dev libboost-all-dev \
        libgl1 libglfw3 libglew2.2 libosmesa6 libegl1 \
    && rm -rf /var/lib/apt/lists/*

# ndcurves (header-only use, no python bindings / tests)
ARG NDCURVES_REF=v1.1.1
RUN git clone --depth 1 --branch ${NDCURVES_REF} --recurse-submodules --shallow-submodules \
        https://github.com/loco-3d/ndcurves.git /tmp/ndcurves \
    && cmake -S /tmp/ndcurves -B /tmp/ndcurves/build -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_PYTHON_INTERFACE=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=/usr/local \
    && cmake --build /tmp/ndcurves/build -j"$(nproc)" && cmake --install /tmp/ndcurves/build \
    && rm -rf /tmp/ndcurves

# Source (acados submodule is fetched below if the checkout doesn't contain it)
ARG ACADOS_REF=04c91bf221b32050050052b9a749e6835f8b64f4
WORKDIR /opt/DWMPC
COPY . .
RUN if [ ! -f third_party/acados/CMakeLists.txt ]; then \
        rm -rf third_party/acados && git clone https://github.com/acados/acados.git third_party/acados \
        && git -C third_party/acados checkout ${ACADOS_REF} \
        && git -C third_party/acados submodule update --init --recursive; \
    fi

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build

ENV ACADOS_SOURCE_DIR=/opt/DWMPC/third_party/acados \
    LD_LIBRARY_PATH=/usr/local/lib:/opt/DWMPC/third_party/acados/lib:/usr/lib/dls2/controllers/dwmpc \
    PYTHONPATH=/usr/lib/dls2/controllers/dwmpc

# Python example deps
RUN git clone --depth 1 https://github.com/iit-DLSLab/gym-quadruped.git /opt/gym-quadruped \
    && pip3 install --no-cache-dir /opt/gym-quadruped

CMD ["python3", "example/sim_dwmpc.py"]
