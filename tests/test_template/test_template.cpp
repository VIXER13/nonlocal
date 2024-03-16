//Здесь добавляем все необходимые библиотеки для тестирования
#include <vector>

#include <boost/ut.hpp>

namespace {
    using namespace boost::ut;
    constexpr double pi = 3.14;
    //здесь пишем весь функционал, который хотим использовать для тестирования

    std::vector<double> template_test(){
        std::vector<double> result;
        constexpr double q = pi / 180.0;
        for (std::size_t i = 0; i <= 180; ++i)
            result.push_back(i * q);
        return result;
    }

    const suite<"template"> _ = [] {
        "vector"_test = [] {
            std::vector<double> res = template_test();
            expect(res.size() > 0);
            for (auto& x: res)
                expect(x == x);
        };
    };
    
}
