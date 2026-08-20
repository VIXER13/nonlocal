#include <boost/ut.hpp>

namespace boost::ut {
    class safe_runner : public runner<reporter_junit<printer>> {
    public:
        [[nodiscard]] auto run(run_cfg rc = {}) -> bool {
            run_ = true;
            reporter_.on(events::run_begin{.argc = rc.argc, .argv = rc.argv});
            for (const auto &[suite, suite_name] : suites_) {
                if constexpr (requires { reporter_.on(events::suite_begin{}); }) {
                    reporter_.on(events::suite_begin{.type = "suite", .name = suite_name});
                }
                constexpr auto type = "placeholder";
                std::string name = std::string(suite_name) + " suite preliminaries";
                try {
                    suite();
                } catch (const std::exception& exception) {
                    ++fails_;
                    reporter_.on(events::test_begin{.type = type, .name = name});
                    reporter_.on(events::exception{exception.what()});
                    reporter_.on(events::test_end{.type = type, .name = name});
                } catch (...) {
                    ++fails_;
                    reporter_.on(events::test_begin{.type = type, .name = name});
                    reporter_.on(events::exception{"Unknown exception"});
                    reporter_.on(events::test_end{.type = type, .name = name});
                }
                if constexpr (requires { reporter_.on(events::suite_end{}); }) {
                    reporter_.on(events::suite_end{.type = "suite", .name = suite_name});
                }
            }
            suites_.clear();

            if (rc.report_errors) {
                report_summary();
            }

            return fails_ > 0;
        }
    };

    template <>
    inline auto cfg<override> = safe_runner{};
}

int main(int argc, const char **argv) {
    return boost::ut::cfg<>.run({.argc = argc, .argv = argv});
}