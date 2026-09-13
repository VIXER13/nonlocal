#include <iostream>
#include <unordered_set>

#include "cli.hpp"

cli::iarg::iarg(cli& cli) : _cli(cli) {}

cli::pos_arg::pos_arg(cli& cli, std::string_view descr, std::size_t pos) : iarg(cli), _pos(pos) {
  _describer = [this, descr](std::ostream& os) {
    os << "arg" << _pos + 1 << "\t<" << _type_name << ">  " << descr;
    if (!_default_repr.empty()) os << " (default: " << _default_repr << ")";
    os << '\n';
  };
}

cli::kv_arg::kv_arg(cli& cli, std::string_view descr, std::string_view name) : iarg(cli), _name(name) {
  _cli.w_name = std::max(_cli.w_name, _name.size());
  _describer = [this, descr](std::ostream& os) {
    os << std::setw(_cli.w_name) << _name << (_alias.empty() ? "  " : ", ") << std::setw(_cli.w_alias)
       << (_alias.empty() ? "" : _alias) << std::setw(_cli.w_type - _type_name.size() + 2) << "<" << _type_name << ">  "
       << descr;
    if (!_default_repr.empty()) os << " (default: " << _default_repr << ")";
    os << '\n';
  };
}

cli::kv_arg& cli::kv_arg::alias(std::string_view name) {
  _cli._named[name] = _cli._named[_name];
  _alias = name;
  _cli.w_alias = std::max(_cli.w_alias, _alias.size());
  return *this;
}

// flag_arg implementations
cli::flag_arg::flag_arg(cli& cli, std::string_view descr, std::string_view name) : kv_arg(cli, descr, name) {
  _parser = [this](std::string_view) { _value = true; };
  _type_name = "flag";
  _value = false;
}

// Cli method implementations
cli::flag_arg& cli::flag(std::string_view name, std::string_view descr) {
  auto arg_ptr = std::make_shared<flag_arg>(*this, descr, name);
  _named[name] = arg_ptr;
  return static_cast<flag_arg&>(*arg_ptr);
}

cli::kv_arg& cli::arg(std::string_view name, std::string_view descr) {
  auto arg_ptr = std::make_shared<kv_arg>(*this, descr, name);
  _named[name] = arg_ptr;
  return *arg_ptr;
}

cli::pos_arg& cli::pos(std::string_view descr) {
  auto arg_ptr = std::make_shared<pos_arg>(*this, descr, _positional.size());
  _positional.emplace_back(arg_ptr);
  return *arg_ptr;
}

void cli::parse(int argc, char* argv[]) {
  auto it = _positional.begin();
  for (int i = 1; i < argc; ++i) {
    std::string_view arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      flag("--help", "Show this help").alias("-h");
      help(argv[0]);
      std::exit(0);
    } else if (arg == "--") {
      while (it != _positional.end())
        (*it++)->_parser(arg);
      i = argc;
    } else if (arg.starts_with("-")) {
      std::string_view key = arg, val = "";
      if (size_t eq = arg.find('='); eq != std::string_view::npos) {
        key = arg.substr(0, eq);
        val = arg.substr(eq + 1);
      }
      _named.at(key)->_parser(val);
    } else {
      (*it++)->_parser(arg);
    }
  }
}

void cli::help(std::string_view executable) const {
  std::cout << "Usage: " << executable;
  if (!_named.empty()) std::cout << " [--option-name=<value>]...";
  std::cout << " [--flag-name]...";
  if (!_positional.empty()) std::cout << " arg1...";
  std::cout << '\n';

  if (!_positional.empty()) {
    std::cout << "\nArgs:\n";
    for (auto&& pos : _positional)
      pos->_describer(std::cout);
  }
  std::cout << "\nOptions:\n";
  std::unordered_set<iarg*> used;
  for (auto&& [name, arg] : _named)
    if (!used.count(arg.get())) {
      arg->_describer(std::cout);
      used.insert(arg.get());
    }
}