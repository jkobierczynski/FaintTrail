// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Process memory use, for profiling.
#include "common.hpp"

#include <cstdio>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/resource.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#endif
#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace fainttrail::detail {

double current_rss_mb() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) return pmc.WorkingSetSize / 1048576.0;
    return 0.0;
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
        return info.resident_size / 1048576.0;
    return 0.0;
#else
    long pages = 0, resident = 0;
    std::FILE* f = std::fopen("/proc/self/statm", "r");
    if (!f) return 0.0;
    const int n = std::fscanf(f, "%ld %ld", &pages, &resident);
    std::fclose(f);
    if (n != 2) return 0.0;
    return static_cast<double>(resident) * static_cast<double>(sysconf(_SC_PAGESIZE)) / 1048576.0;
#endif
}

double peak_rss_mb() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) return pmc.PeakWorkingSetSize / 1048576.0;
    return 0.0;
#else
    rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
#if defined(__APPLE__)
    return ru.ru_maxrss / 1048576.0;  // bytes
#else
    return ru.ru_maxrss / 1024.0;     // kilobytes
#endif
#endif
}

void release_free_memory() {
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
}

}  // namespace fainttrail::detail

namespace fainttrail {

void tune_allocator_for_images() {
#if defined(__GLIBC__)
    // glibc raises its mmap threshold up to 32 MB after large blocks are freed, so
    // that later frame-sized buffers come from the heap and freed ones are kept by
    // the process. A fixed 1 MB threshold returns every frame buffer on release.
    mallopt(M_MMAP_THRESHOLD, 1 << 20);
#endif
}

}  // namespace fainttrail
