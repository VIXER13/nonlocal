#include <array>
#include <exception>
#include <stacktrace>

#include "runner.hpp"

// Storage to store up to 5 nested stack traces for exceptions
//   should be enough for most cases but can be increased if needed
thread_local std::array<std::stacktrace, 5> s_stacktraces;

// Forward declaration for the original __cxa_throw function
//   that will be called from our wrapper
extern "C" auto __real___cxa_throw(void* thrown_object, std::type_info* tinfo,
                                   void (*dest)(void*)) -> void;

// Wrapper function for __cxa_throw that will be called
//   instead of the original one
extern "C" auto __wrap___cxa_throw(void* thrown_object, std::type_info* tinfo,
                                   void (*dest)(void*)) -> void {
  // Called when an exception is thrown, at the call site of the `throw`

  // std::uncaught_exceptions() returns the number of currently
  //   active exceptions that have been thrown but not yet caught
  auto exception_count = std::uncaught_exceptions();

  // If there's still some space on the exception stack,
  //   capture a new stacktrace and add it to the stack
  if (exception_count < ssize_t(s_stacktraces.size())) {
    // Add the current stack trace to the end of the
    //   stack trace storage array skipping the current
    //   frame which is the __wrap___cxa_throw function itself
    s_stacktraces[exception_count] = std::stacktrace::current(1);
  }

  // Forward to original __cxa_throw()
  __real___cxa_throw(thrown_object, tinfo, dest);
}

// Simple helper function to print the stack trace
//   of the exception that's currently being caught
extern "C" auto get_stacktrace() -> std::string {
  auto exception_count = std::uncaught_exceptions();
  if (exception_count < 0) {
    return "No active exception";
  } else if (exception_count >= s_stacktraces.size()) {
    return "Too many active exceptions";
  }
  return std::to_string(s_stacktraces[exception_count]);
}

int main(int argc, const char** argv) {
  return boost::ut::cfg<>.run({.argc = argc, .argv = argv});
}