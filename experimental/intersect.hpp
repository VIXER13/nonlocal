// Вычисление ориентированной площади треугольника на точках a, b, c.
template<class T>
T oriented_triangle_area(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c) noexcept {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

// Проверка отрезков AB и CD на предмет пересечения.
template<class T>
bool intersect(const std::array<T, 2>& a, const std::array<T, 2>& b, const std::array<T, 2>& c, const std::array<T, 2>& d) noexcept {
    static constexpr auto check_segment = [](T a, T b, T c, T d) {
        if (a > b) std::swap(a, b);
        if (c > d) std::swap(c, d);
        return std::max(a, c) <= std::min(b, d);
    };
    
    return check_segment(a[0], b[0], c[0], d[0]) && 
           check_segment(a[1], b[1], c[1], d[1]) &&
           oriented_triangle_area(a, b, c) * oriented_triangle_area(a, b, d) <= 0 &&
           oriented_triangle_area(c, d, a) * oriented_triangle_area(c, d, b) <= 0;
}

template<class T, class I>
bool mesh_2d<T, I>::check_intersect(const std::array<T, 2>& A, const std::array<T, 2>& B) const {
    for(size_t b = 0; b < boundary_groups_count(); ++b) {
        for(size_t el = 0; el < elements_count(b); ++el) {
            const auto& be = element_1d(element_1d_type(b, el));
            for(size_t i = 0; i < be->nodes_count()-1; ++i) {
                if(intersect(node(node_number(b, el, i)), node(node_number(b, el, i+1)), A, B))
                    return true;
            }
        }
    }
    return false;
}