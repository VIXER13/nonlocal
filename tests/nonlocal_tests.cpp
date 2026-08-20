#include <boost/ut.hpp>

#include <csignal>
#include <cstdlib>

namespace {
    void signal_handler(const int sig) {
        std::signal(sig, SIG_DFL);
        std::string err = std::string(boost::ut::colors{}.fail);
        err += "FAILED\nSignal received: ";
        err += [sig] {
            switch(sig) {
            case SIGABRT: return "SIGABRT";
            case SIGSEGV: return "SIGSEGV";
            case SIGFPE:  return "SIGFPE";
            case SIGILL:  return "SIGILL";
            case SIGTERM: return "SIGTERM";
#ifndef _WIN32
            case SIGBUS:  return "SIGBUS";
#endif
            default:      return "Unknown signal";
            }   
        }();
        err += boost::ut::colors{}.none;
        boost::ut::cfg<>.on(boost::ut::events::log{err});
        boost::ut::cfg<>.on(boost::ut::events::fatal_assertion{});
    }

    void install_signal_handlers() {
        std::signal(SIGABRT, signal_handler);
        std::signal(SIGSEGV, signal_handler);
        std::signal(SIGFPE, signal_handler);
        std::signal(SIGILL, signal_handler);
        std::signal(SIGTERM, signal_handler);
#ifndef _WIN32
        std::signal(SIGBUS, signal_handler);
#endif
    }
}

int main() {
    install_signal_handlers();
}