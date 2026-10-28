// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2026 Jurgen Kobierczynski
//
// Starting external programs (ASTAP, curl) without a shell.
#include "process.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace fainttrail::detail {

// Run a program with arguments, wait up to timeout_s. Returns the exit code, -1
// if it could not be started, -2 on timeout.
int run_program(const std::string& exe, const std::vector<std::string>& args, int timeout_s, const std::string& log_path) {
#ifdef _WIN32
    auto quote = [](const std::string& a) {
        if (!a.empty() && a.find_first_of(" \t\"") == std::string::npos) return a;
        std::string q = "\"";
        for (char c : a) {
            if (c == '"') q += '\\';
            q += c;
        }
        return q + "\"";
    };
    std::string cmd = quote(exe);
    for (const auto& a : args) cmd += " " + quote(a);
    auto widen = [](const std::string& s) {
        const int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, nullptr, 0);
        std::wstring w(static_cast<std::size_t>(n > 0 ? n : 1), L'\0');
        if (n > 0) MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, &w[0], n);
        return w;
    };
    std::wstring wcmd = widen(cmd), wexe = widen(exe), wlog = widen(log_path);
    SECURITY_ATTRIBUTES sa{sizeof(sa), nullptr, TRUE};
    HANDLE logh = CreateFileW(wlog.c_str(), GENERIC_WRITE, FILE_SHARE_READ, &sa, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    STARTUPINFOW si{};
    si.cb = sizeof(si);
    if (logh != INVALID_HANDLE_VALUE) {
        si.dwFlags = STARTF_USESTDHANDLES;
        si.hStdOutput = logh;
        si.hStdError = logh;
        si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    }
    PROCESS_INFORMATION pi{};
    if (!CreateProcessW(wexe.c_str(), &wcmd[0], nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi)) {
        if (logh != INVALID_HANDLE_VALUE) CloseHandle(logh);
        return -1;
    }
    const DWORD w = WaitForSingleObject(pi.hProcess, static_cast<DWORD>(timeout_s) * 1000u);
    int code = -2;
    if (w == WAIT_OBJECT_0) {
        DWORD c = 0;
        GetExitCodeProcess(pi.hProcess, &c);
        code = static_cast<int>(c);
    } else {
        TerminateProcess(pi.hProcess, 1);
    }
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    if (logh != INVALID_HANDLE_VALUE) CloseHandle(logh);
    return code;
#else
    const pid_t pid = fork();
    if (pid < 0) return -1;
    if (pid == 0) {
        if (std::FILE* f = std::fopen(log_path.c_str(), "w")) {
            dup2(fileno(f), 1);
            dup2(fileno(f), 2);
        }
        std::vector<char*> argv;
        argv.push_back(const_cast<char*>(exe.c_str()));
        for (const auto& a : args) argv.push_back(const_cast<char*>(a.c_str()));
        argv.push_back(nullptr);
        execv(exe.c_str(), argv.data());
        _exit(127);
    }
    const auto t0 = std::chrono::steady_clock::now();
    for (;;) {
        int st = 0;
        const pid_t r = waitpid(pid, &st, WNOHANG);
        if (r == pid) {
            if (WIFEXITED(st)) return WEXITSTATUS(st) == 127 ? -1 : WEXITSTATUS(st);
            return -3;
        }
        if (r < 0) return -1;
        if (std::chrono::steady_clock::now() - t0 > std::chrono::seconds(timeout_s)) {
            kill(pid, SIGKILL);
            waitpid(pid, &st, 0);
            return -2;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
#endif
}


std::string find_program(const std::vector<std::string>& names, const std::vector<std::string>& extra_dirs) {
#ifdef _WIN32
    const char sep = ';';
#else
    const char sep = ':';
#endif
    std::vector<std::string> dirs = extra_dirs;
    if (const char* path = std::getenv("PATH")) {
        std::string p = path;
        std::size_t s = 0;
        while (s <= p.size()) {
            const std::size_t e = p.find(sep, s);
            const std::string d = p.substr(s, e == std::string::npos ? std::string::npos : e - s);
            if (!d.empty()) dirs.push_back(d);
            if (e == std::string::npos) break;
            s = e + 1;
        }
    }
    std::error_code ec;
    for (const auto& n : names)
        for (const auto& d : dirs) {
            const fs::path c = fs::path(d) / n;
            if (fs::is_regular_file(c, ec)) return c.string();
        }
    return {};
}

}  // namespace fainttrail::detail
