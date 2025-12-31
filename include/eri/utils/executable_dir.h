#pragma once

#include <unistd.h>
#include <limits.h>

#include <string>
#include <stdexcept>

namespace eri::utils {

// Resolve directory of the running executable (Linux)
static std::string executable_dir() {
    char buffer[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);

    if (len == -1)
        throw std::runtime_error("Failed to resolve /proc/self/exe");

    buffer[len] = '\0';

    std::string path(buffer);
    return path.substr(0, path.find_last_of('/'));
}

}