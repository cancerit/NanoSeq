# Builder stage
FROM ubuntu:24.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    software-properties-common \
    ca-certificates \
    lsb-release \
    && rm -rf /var/lib/apt/lists/*
RUN wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc \
    && echo "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-20 main" > /etc/apt/sources.list.d/llvm.list
RUN apt-get update && apt-get install -y \
    clang-20 \
    lldb-20 \
    lld-20 \
    libc++-20-dev \
    libc++abi-20-dev \
    libclang-20-dev \
    libclang-common-20-dev \
    && rm -rf /var/lib/apt/lists/*
RUN update-alternatives --install /usr/bin/clang clang /usr/bin/clang-20 100 \
    && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-20 100 \
    && update-alternatives --install /usr/bin/lld lld /usr/bin/lld-20 100
RUN apt-get update && apt-get install -y \
    libhts-dev \
    libdeflate-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y \
    make \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*
ENV CC=clang-20
ENV CXX=clang++-20
# -fexperimental-library is still needed for <flat_map> in libc++20.
# On Debian/Ubuntu apt packages clang does NOT auto-add /usr/lib/llvm-N/lib
# to the linker search path, so libc++experimental.a is never found.
# The -L flag below fixes that.
ENV CXXFLAGS="-std=c++23 -stdlib=libc++ -fexperimental-library"
ENV LDFLAGS="-stdlib=libc++ -fexperimental-library -L/usr/lib/llvm-20/lib"
WORKDIR /workspace
COPY build.sh ./
COPY src/ ./src
RUN ./build.sh

# Runtime stage
FROM ubuntu:24.04 AS runtime
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    ca-certificates \
    && wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc \
    && echo "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-20 main" > /etc/apt/sources.list.d/llvm.list \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y --no-install-recommends \
    libc++1-20 \
    libc++abi1-20 \
    libhts3 \
    libdeflate0 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libssl3 \
    strace \
    && apt-get purge -y wget gnupg \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
COPY --from=builder /workspace/dsa ./dsa
USER ubuntu
ENTRYPOINT ["/workspace/dsa"]
