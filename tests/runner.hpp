#pragma once
#include <boost/ut.hpp>
#include <typeinfo>

// Forward declaration for the original __cxa_throw function
//   that will be called from our wrapper
extern "C" auto __real___cxa_throw(void* thrown_object, std::type_info* tinfo,
                                   void (*dest)(void*)) -> void;

// Wrapper function for __cxa_throw that will be called
//   instead of the original one
extern "C" auto __wrap___cxa_throw(void* thrown_object, std::type_info* tinfo,
                                   void (*dest)(void*)) -> void;

// Simple helper function to print the stack trace
//   of the exception that's currently being caught
extern "C" auto get_stacktrace() -> std::string;

namespace boost::inline ext::ut::inline v2_3_1 {
template <class TReporter = reporter<printer>, auto MaxPathSize = 16>
class safe_runner {
  class filter {
    static constexpr auto delim = ".";

   public:
    constexpr /*explicit(false)*/ filter(std::string_view _filter = {})
        : path_{utility::split(_filter, delim)} {}

    template <class TPath>
    constexpr auto operator()(const std::size_t level, const TPath& path) const
        -> bool {
      for (auto i = 0u; i < math::min_value(level + 1, std::size(path_)); ++i) {
        if (not utility::is_match(path[i], path_[i])) {
          return false;
        }
      }
      return true;
    }

   private:
    std::vector<std::string_view> path_{};
  };

 public:
  constexpr safe_runner() = default;
  constexpr safe_runner(TReporter reporter, std::size_t suites_size)
      : reporter_{std::move(reporter)}, suites_(suites_size) {}

  ~safe_runner() {
    const auto should_run = not run_;

    if (should_run) {
      static_cast<void>(run());
    }

    if (not dry_run_) {
      report_summary();
    }

    if (should_run and fails_) {
      std::exit(-1);
    }
  }

  auto operator=(const options& options) {
    filter_ = options.filter;
    tag_ = options.tag;
    dry_run_ = options.dry_run;
    reporter_ = {options.colors};
  }

  template <class TSuite>
  auto on(events::suite<TSuite> suite) {
    suites_.emplace_back(suite.run, suite.name);
  }

  template <class... Ts>
  auto on(events::test<Ts...> test) {
    path_[level_] = test.name;

    if (detail::cfg::list_tags) {
      std::for_each(test.tag.cbegin(), test.tag.cend(), [](const auto& tag) {
        std::cout << "tag: " << tag << std::endl;
      });
      return;
    }

    auto execute = std::empty(test.tag);
    for (const auto& tag_element : test.tag) {
      if (utility::is_match(tag_element, "skip") && !detail::cfg::show_tests &&
          !detail::cfg::show_test_names) {
        on(events::skip<>{.type = test.type, .name = test.name});
        return;
      }

      for (const auto& ftag : tag_) {
        if (utility::is_match(tag_element, ftag)) {
          execute = true;
          break;
        }
      }
    }

    if (!detail::cfg::query_pattern.empty()) {
      const static auto regex = detail::cfg::query_regex_pattern;
      bool matches = utility::regex_match(test.name.data(), regex.c_str());
      for (const auto& tag2 : test.tag) {
        matches |= utility::regex_match(tag2.data(), regex.c_str());
      }
      if (matches) {
        execute = !detail::cfg::invert_query_pattern;
      } else {
        execute = detail::cfg::invert_query_pattern;
      }
    }

    if (detail::cfg::show_tests || detail::cfg::show_test_names) {
      if (!detail::cfg::show_test_names) {
        std::cout << "matching test: ";
      }
      std::cout << test.name << std::endl;
      return;
    }

    if (not execute) {
      on(events::skip<>{.type = test.type, .name = test.name});
      return;
    }

    if (filter_(level_, path_)) {
      if (not level_++) {
        reporter_.on(events::test_begin{
            .type = test.type, .name = test.name, .location = test.location});
      } else {
        reporter_.on(events::test_run{.type = test.type, .name = test.name});
      }

      if (dry_run_) {
        for (auto i = 0u; i < level_; ++i) {
          std::cout << (i ? "." : "") << path_[i];
        }
        std::cout << '\n';
      }

#if defined(__cpp_exceptions)
      try {
#endif
        test();
#if defined(__cpp_exceptions)
      } catch (const events::fatal_assertion&) {
      } catch (const std::exception& exception) {
        ++fails_;
        std::string msg = std::string(exception.what()) + "\nStacktrace:\n" +
                          get_stacktrace();
        reporter_.on(events::exception{msg.c_str()});
      } catch (...) {
        ++fails_;
        std::string msg = "Unknown exception\nStacktrace:\n" + get_stacktrace();
        reporter_.on(events::exception{msg.c_str()});
      }
#endif

      if (not --level_) {
        reporter_.on(events::test_end{.type = test.type, .name = test.name});
      } else {  // N.B. prev. only root-level tests were signalled on finish
        if constexpr (requires {
                        reporter_.on(events::test_finish{.type = test.type,
                                                         .name = test.name});
                      }) {
          reporter_.on(
              events::test_finish{.type = test.type, .name = test.name});
        }
      }
    }
  }

  template <class... Ts>
  auto on(events::skip<Ts...> test) {
    reporter_.on(events::test_skip{.type = test.type, .name = test.name});
  }

  template <class TExpr>
  [[nodiscard]] auto on(events::assertion<TExpr> assertion) -> bool {
    if (dry_run_) {
      return true;
    }

    if (static_cast<bool>(assertion.expr)) {
      reporter_.on(events::assertion_pass<TExpr>{
          .expr = assertion.expr, .location = assertion.location});
      return true;
    }

    ++fails_;
    reporter_.on(events::assertion_fail<TExpr>{.expr = assertion.expr,
                                               .location = assertion.location});
    return false;
  }

  auto on(events::fatal_assertion fatal_assertion) {
    reporter_.on(fatal_assertion);

#if defined(__cpp_exceptions)
    if (not level_) {
      report_summary();
    }
    throw fatal_assertion;
#else
    if (level_) {
      reporter_.on(events::test_end{});
    }
    report_summary();
    std::abort();
#endif
  }

  template <class TMsg>
  auto on(events::log<TMsg> l) {
    reporter_.on(l);
  }

  [[nodiscard]] auto run(run_cfg rc = {}) -> bool {
    run_ = true;
    reporter_.on(events::run_begin{.argc = rc.argc, .argv = rc.argv});
    for (const auto& [suite, suite_name] : suites_) {
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

  auto report_summary() -> void {
    if (static auto once = true; once) {
      once = false;
      reporter_.on(events::summary{});
    }
  }

 protected:
  TReporter reporter_{};
  std::vector<std::pair<void (*)(), std::string_view>> suites_{};
  std::size_t level_{};
  bool run_{};
  std::size_t fails_{};
  std::array<std::string_view, MaxPathSize> path_{};
  filter filter_{};
  std::vector<std::string_view> tag_{};
  bool dry_run_{};
};

template <>
inline auto cfg<override> = safe_runner<reporter_junit<printer>>{};
}