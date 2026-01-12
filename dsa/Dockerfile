# Builder stage
FROM ubuntu:24.04 AS builder

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Update and install basic dependencies
RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    software-properties-common \
    ca-certificates \
    lsb-release \
    && rm -rf /var/lib/apt/lists/*

# Add LLVM apt repository
RUN wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc \
    && echo "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-18 main" > /etc/apt/sources.list.d/llvm.list

# Install LLVM 18 with all components including libc++
RUN apt-get update && apt-get install -y \
    clang-18 \
    lldb-18 \
    lld-18 \
    libc++-18-dev \
    libc++abi-18-dev \
    libclang-18-dev \
    libclang-common-18-dev \
    && rm -rf /var/lib/apt/lists/*

# Create symlinks for convenience
RUN update-alternatives --install /usr/bin/clang clang /usr/bin/clang-18 100 \
    && update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-18 100 \
    && update-alternatives --install /usr/bin/lld lld /usr/bin/lld-18 100

# Install htslib dependencies
RUN apt-get update && apt-get install -y \
    libhts-dev \
    libdeflate-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libpthread-stubs0-dev \
    && rm -rf /var/lib/apt/lists/*

# Install build tools
RUN apt-get update && apt-get install -y \
    make \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for compilation
ENV CC=clang-18
ENV CXX=clang++-18
ENV CXXFLAGS="-stdlib=libc++ -fexperimental-library"
ENV LDFLAGS="-stdlib=libc++ -fexperimental-library"

# Set working directory
WORKDIR /workspace

# Copy source and build
COPY build.sh ./
COPY src/ ./src
RUN ./build.sh

# Runtime stage
FROM ubuntu:24.04 AS runtime

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Add LLVM apt repository (needed for libc++ runtime)
RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    ca-certificates \
    && wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc \
    && echo "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-18 main" > /etc/apt/sources.list.d/llvm.list \
    && rm -rf /var/lib/apt/lists/*

# Install only runtime libraries (no -dev packages)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libc++1-18 \
    libc++abi1-18 \
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

# Set working directory
WORKDIR /workspace

# Copy the compiled binary from the builder stage
COPY --from=builder /workspace/dsa ./dsa

# Switch to non-root user
USER ubuntu

ENTRYPOINT ["/workspace/dsa"]
