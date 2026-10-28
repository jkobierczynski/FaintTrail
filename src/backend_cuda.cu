// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// CUDA stacking backend. One thread per output pixel; the per-pixel arithmetic
// is cuda_core::stack_pixel, shared with the host emulation the unit tests run.
// Frames are summed in ascending order in single precision and the result is
// bit-identical to the CPU backend (built with --fmad=false; the kernel has no
// multiply-add to contract anyway, the flag guards against future changes).
//
// Flow per pass: plan (host) -> for each row band: upload the frame rows the band
// needs -> kernel launches over chunks of motions, each thread carrying its
// running maximum in device memory between launches -> download. Progress and
// cancellation are checked between launches.
#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

#include "cuda_plan.hpp"
#include "cuda_stack_core.hpp"

namespace fainttrail::detail {

namespace {

void check(cudaError_t e, const char* what) {
    if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA: ") + what + ": " + cudaGetErrorString(e));
}

template <typename T>
struct DeviceBuffer {
    T* p = nullptr;
    std::size_t n = 0;
    DeviceBuffer() = default;
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    ~DeviceBuffer() {
        if (p) cudaFree(p);
    }
    void reserve(std::size_t count, const char* what) {
        if (count <= n) return;
        if (p) cudaFree(p);
        p = nullptr;
        n = 0;
        check(cudaMalloc(reinterpret_cast<void**>(&p), std::max<std::size_t>(1, count) * sizeof(T)), what);
        n = count;
    }
};

__global__ void stack_kernel(cuda_core::BandView v, int y0, int y1, int j0, int j1, float* best, int* arg) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = y0 + blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= v.w || y >= y1) return;
    const std::size_t o = static_cast<std::size_t>(y - y0) * v.w + x;
    float b;
    int a;
    if (j0 == 0) {
        b = __int_as_float(0xff800000);  // -inf (the INFINITY macro of MSVC overflows a float)
        a = -1;
    } else {
        b = best[o];
        a = arg[o];
    }
    cuda_core::stack_pixel(v, x, y, j0, j1, b, a);
    best[o] = b;
    arg[o] = a;
}

int resolve_device(int requested) {
    if (requested >= 0) return requested;
    if (const char* e = std::getenv("FAINTTRAIL_CUDA_DEVICE")) return std::atoi(e);
    return 0;
}

class CudaBackend final : public SearchBackend {
public:
    explicit CudaBackend(const CudaOptions& o) : opt_(o), device_(resolve_device(o.device)) {
        int count = 0;
        check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
        if (device_ < 0 || device_ >= count)
            throw std::runtime_error("CUDA device " + std::to_string(device_) + " not present (" + std::to_string(count) + " found)");
        check(cudaGetDeviceProperties(&prop_, device_), "cudaGetDeviceProperties");
    }

    std::string name() const override { return "cuda"; }

    std::string description() const override {
        char buf[256];
        std::snprintf(buf, sizeof buf, "CUDA device %d: %s, %.1f GB, compute capability %d.%d, %d SMs", device_, prop_.name,
                      prop_.totalGlobalMem / 1073741824.0, prop_.major, prop_.minor, prop_.multiProcessorCount);
        return buf;
    }

    StackResult run(const StackPass& pass, const ProgressFn& progress) override {
        using clk = std::chrono::steady_clock;
        const auto t_start = clk::now();
        check(cudaSetDevice(device_), "cudaSetDevice");
        const auto& frames = *pass.frames;
        validate_pass(pass);
        const int w = frames[0].width, h = frames[0].height, n = static_cast<int>(frames.size());

        // device memory: everything except the frames is small and per band
        std::size_t free_b = 0, total_b = 0;
        check(cudaMemGetInfo(&free_b, &total_b), "cudaMemGetInfo");
        const double nv_d = static_cast<double>(pass_motions(pass));
        const double fixed = 8.0 * nv_d * n + 4.0 * (n + 1) + 64.0 * 1048576.0;  // shifts, 1/sqrt, slack
        double budget = 0.8 * static_cast<double>(free_b) - fixed;
        if (opt_.max_mb > 0) budget = std::min(budget, opt_.max_mb * 1048576.0);
        // the output of a band (best + arg per pixel) shares the budget with its frame
        // rows; per output row it costs 8 w bytes against 4 w n for the frames
        const double frame_share = static_cast<double>(n) / (n + 2.0);
        const CudaPlan p = make_cuda_plan(pass, std::max(1.0, budget * frame_share), opt_.launch_adds);

        DeviceBuffer<int> d_dx, d_dy;
        DeviceBuffer<float> d_inv, d_frames, d_best;
        DeviceBuffer<int> d_arg;
        d_dx.reserve(p.dx.size(), "allocate shifts");
        d_dy.reserve(p.dy.size(), "allocate shifts");
        d_inv.reserve(p.inv_sqrt.size(), "allocate 1/sqrt table");
        check(cudaMemcpy(d_dx.p, p.dx.data(), p.dx.size() * sizeof(int), cudaMemcpyHostToDevice), "upload shifts");
        check(cudaMemcpy(d_dy.p, p.dy.data(), p.dy.size() * sizeof(int), cudaMemcpyHostToDevice), "upload shifts");
        check(cudaMemcpy(d_inv.p, p.inv_sqrt.data(), p.inv_sqrt.size() * sizeof(float), cudaMemcpyHostToDevice), "upload table");

        StackResult out;
        out.max_snr = Image(w, h, 0.0f);
        out.best_v.assign(static_cast<std::size_t>(w) * h, -1);
        std::vector<float> best;
        std::vector<int> arg;
        double t_up = 0, t_kernel = 0, t_down = 0;
        const double total_work = static_cast<double>(p.bands.size()) * p.nv;
        double work_done = 0;
        const dim3 block(32, 8);

        for (const auto& b : p.bands) {
            const int rows = b.r1 - b.r0, out_rows = b.y1 - b.y0;
            const std::size_t np = static_cast<std::size_t>(out_rows) * w;
            auto t0 = clk::now();
            d_frames.reserve(static_cast<std::size_t>(n) * rows * w, "allocate frames");
            d_best.reserve(np, "allocate output");
            d_arg.reserve(np, "allocate output");
            for (int i = 0; i < n; ++i)
                check(cudaMemcpy(d_frames.p + static_cast<std::size_t>(i) * rows * w,
                                 frames[static_cast<std::size_t>(i)].data.data() + static_cast<std::size_t>(b.r0) * w,
                                 static_cast<std::size_t>(rows) * w * sizeof(float), cudaMemcpyHostToDevice),
                      "upload frames");
            auto t1 = clk::now();
            t_up += std::chrono::duration<double>(t1 - t0).count();

            cuda_core::BandView v;
            v.frames = d_frames.p;
            v.w = w;
            v.h = h;
            v.r0 = b.r0;
            v.rows = rows;
            v.n = n;
            v.dx = d_dx.p;
            v.dy = d_dy.p;
            v.inv_sqrt = d_inv.p;
            v.min_count = p.min_count;
            const dim3 grid((w + block.x - 1) / block.x, (out_rows + block.y - 1) / block.y);
            for (int j0 = 0; j0 < p.nv; j0 += p.chunk) {
                const int j1 = std::min(p.nv, j0 + p.chunk);
                stack_kernel<<<grid, block>>>(v, b.y0, b.y1, j0, j1, d_best.p, d_arg.p);
                check(cudaGetLastError(), "launch stack kernel");
                check(cudaDeviceSynchronize(), "stack kernel");
                work_done += j1 - j0;
                if (progress && !progress(work_done / total_work, "stack")) throw Cancelled();
            }
            auto t2 = clk::now();
            t_kernel += std::chrono::duration<double>(t2 - t1).count();

            best.resize(np);
            arg.resize(np);
            check(cudaMemcpy(best.data(), d_best.p, np * sizeof(float), cudaMemcpyDeviceToHost), "download");
            check(cudaMemcpy(arg.data(), d_arg.p, np * sizeof(int), cudaMemcpyDeviceToHost), "download");
            for (std::size_t o = 0; o < np; ++o)
                if (arg[o] >= 0) {
                    out.max_snr.data[static_cast<std::size_t>(b.y0) * w + o] = best[o];
                    out.best_v[static_cast<std::size_t>(b.y0) * w + o] = arg[o];
                }
            t_down += std::chrono::duration<double>(clk::now() - t2).count();
        }
        const double total = std::chrono::duration<double>(clk::now() - t_start).count();
        out.timings.push_back({"upload", t_up});
        out.timings.push_back({"kernel", t_kernel});
        out.timings.push_back({"download", t_down});
        out.timings.push_back({"setup", std::max(0.0, total - t_up - t_kernel - t_down)});
        if (progress) progress(1.0, "stack");
        return out;
    }

private:
    CudaOptions opt_;
    int device_ = 0;
    cudaDeviceProp prop_{};
};

}  // namespace

int cuda_device_count(std::string* why) {
    int count = 0;
    const cudaError_t e = cudaGetDeviceCount(&count);
    if (e != cudaSuccess) {
        if (why) *why = cudaGetErrorString(e);
        cudaGetLastError();  // clear the sticky error
        return 0;
    }
    if (count == 0 && why) *why = "no CUDA device found";
    return count;
}

std::unique_ptr<SearchBackend> make_cuda_backend(const CudaOptions& options) {
    return std::make_unique<CudaBackend>(options);
}

}  // namespace fainttrail::detail
