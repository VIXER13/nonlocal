#pragma once

#include <any>
#include <charconv>
#include <format>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

template<class T>
constexpr std::string_view type_name() {
    auto loc = std::source_location::current();
    std::string_view func{ loc.function_name() };
    constexpr std::string_view pattern = "[with T = ";
    func.remove_prefix(func.find(pattern) + pattern.size());
    return func.substr(0, func.find_first_of(";]"));
}
template<>
constexpr std::string_view type_name<std::string>() {
    return "string";
}

class cli {
    template<typename T>
    class arg_proxy;

    class iarg {
      protected:
        cli& _cli;
        std::string _type_name, _default_repr;
        std::function<void(std::ostream&)> _describer;
        std::function<void(std::string_view)> _parser;
        std::any _value;

      public:
        iarg(cli& cli);
        virtual ~iarg() = default;

        template<typename T>
        void set_default_parser() {
            if constexpr (std::is_same_v<T, std::string>) {
                _parser = [this](std::string_view repr) { _value = std::string(repr); };
            } else if constexpr (std::is_same_v<T, bool>) {
                _parser = [this](std::string_view repr) {
                    if (repr == "true" || repr == "1" || repr == "yes")
                        _value = true;
                    else if (repr == "false" || repr == "0" || repr == "no")
                        _value = false;
                    else
                        throw std::runtime_error("Unsupported boolean repr");
                };
            } else if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
                _parser = [this](std::string_view repr) {
                    T val{};
                    auto [ptr, ec] = std::from_chars(repr.data(), repr.data() + repr.size(), val);
                    if (ec == std::errc{})
                        _value = val;
                    else
                        throw std::runtime_error("Invalid number");
                };
            }
        }

        template<typename F>
        iarg& cast(F&& func) {
            _parser = [this, f = std::forward<F>(func)](std::string_view repr) { _value = f(repr); };
            return *this;
        }
        template<typename T>
        iarg& dflt(const T& value) {
            if constexpr (std::is_convertible_v<T, std::string>)
                _value = std::string(value);
            else
                _value = value;
            _default_repr = std::format("{}", value);
            return *this;
        }
        template<typename T>
        arg_proxy<T> type() {
            _type_name = type_name<T>();
            _cli.w_type = std::max(_cli.w_type, _type_name.size());
            if (!_parser)
                set_default_parser<T>();
            return arg_proxy<T>(*this);
        }

      protected:
        friend class cli;
        template<typename T>
        friend class arg_proxy;
    };

    template<typename T>
    struct arg_proxy {
        const iarg& _arg;
        const T& operator*() {
            return std::any_cast<const T&>(_arg._value);
        }
        bool has_value() {
            return _arg._value.has_value();
        }
    };

    class pos_arg : public iarg {
        std::size_t _pos;

      public:
        pos_arg(cli& cli, std::string_view descr, std::size_t pos);
    };

    class kv_arg : public iarg {
      protected:
        std::string_view _name;
        std::string_view _alias;

      public:
        kv_arg(cli& cli, std::string_view descr, std::string_view name);
        kv_arg& alias(std::string_view name);
    };

    class flag_arg : public kv_arg {
      public:
        flag_arg(cli& cli, std::string_view descr, std::string_view name);
    };

    std::vector<std::shared_ptr<iarg>> _positional;
    std::unordered_map<std::string_view, std::shared_ptr<iarg>> _named;
    std::size_t w_name = 0, w_alias = 0, w_type = 0;

  public:
    void parse(int argc, char* argv[]);
    void help(std::string_view executable) const;

    kv_arg& arg(std::string_view name, std::string_view descr);
    pos_arg& pos(std::string_view descr);
    flag_arg& flag(std::string_view name, std::string_view descr);
};
