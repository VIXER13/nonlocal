# Stacktrace support library

include(CheckLinkerFlag)
add_library(stacktrace_support INTERFACE)
check_linker_flag(CXX "-lstdc++exp" HAVE_STDCXXEXP)
check_linker_flag(CXX "-lstdc++_libbacktrace" HAVE_STDCXX_LIBBACKTRACE)
if(HAVE_STDCXXEXP)
    target_link_libraries(stacktrace_support INTERFACE stdc++exp)
elseif(HAVE_STDCXX_LIBBACKTRACE)
    target_link_libraries(stacktrace_support INTERFACE stdc++_libbacktrace)
else()
    message(WARNING "Stacktrace: neither stdc++exp nor stdc++_libbacktrace found.")
endif()