#pragma once

#include <array>
#include <concepts>
#include <memory>
#include <type_traits>
#include <vector>

namespace metamath::types {

template<class T>
struct is_array : std::is_array<T> {};

template<class T, size_t N>
struct is_array<std::array<T, N>> : std::true_type {};

template<class T>
constexpr bool is_array_v = is_array<T>::value;

static_assert(!is_array_v<int>, "int is not an array");
static_assert(is_array_v<int[5]>, "int[5] is an array");
static_assert(is_array_v<std::array<int, 5>>, "std::array<int, 5> is an array");

template<class T>
struct container_type : std::type_identity<T> {};

template<class T, size_t N>
struct container_type<std::array<T, N>> : std::type_identity<T> {};

template<class T>
struct container_type<std::vector<T>> : std::type_identity<T> {};

template<class T>
using container_type_t = typename container_type<T>::type;

static_assert(std::is_same_v<container_type_t<int>, int>, "container_type_t<int> should be int");
static_assert(std::is_same_v<container_type_t<std::array<int, 5>>, int>, "container_type_t<std::array<int, 5>> should be int");
static_assert(std::is_same_v<container_type_t<std::vector<int>>, int>, "container_type_t<std::vector<int>> should be int");

template<class T> 
concept arithmetic = std::integral<T> || std::floating_point<T>;

template<class T>
concept copyable = requires(const T& v) { { v.copy() } -> std::convertible_to<std::unique_ptr<T>>; };

}