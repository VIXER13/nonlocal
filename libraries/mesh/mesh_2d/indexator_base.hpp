#pragma once

namespace nonlocal::mesh {

class indexator_base {
    const bool _is_symmetric = false;

protected:
    explicit indexator_base(const bool is_symmetric) noexcept
        : _is_symmetric{is_symmetric} {}
    virtual ~indexator_base() noexcept = default;

public:
    bool is_symmetric() const noexcept {
        return _is_symmetric;
    }

    virtual void reset(const size_t node) = 0;
};

}