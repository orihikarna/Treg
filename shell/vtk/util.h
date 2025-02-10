#pragma once

#include <sys/resource.h>
#include <sys/time.h>

#include <iostream>
#include <string>

long dumpMemoryUsage(const std::string &title = "") {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    const int mem_MB = std::round(usage.ru_maxrss / 1024.0 / 1024.0);
    std::cout << "===[ Memory: " << mem_MB << " MB (" << title << ") ]===" << std::endl;
    return usage.ru_maxrss;  // bytes
  } else
    return 0;
}
