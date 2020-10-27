// Производные высокого порядка для кубических серендиповых элементов

/*
    // Производные высших порядков экспериментальны и не подвергались оптимизации и тестированию.
    Type Nxi2   (const size_t i, const Type xi, const Type eta) const { return basic_Nxi2   [i](xi, eta); }
    Type Nxieta (const size_t i, const Type xi, const Type eta) const { return basic_Nxieta [i](xi, eta); }
    Type Neta2  (const size_t i, const Type xi, const Type eta) const { return basic_Neta2  [i](xi, eta); }
    Type Nxi3   (const size_t i, const Type xi, const Type eta) const { return basic_Nxi3   [i](xi, eta); }
    Type Nxi2eta(const size_t i, const Type xi, const Type eta) const { return basic_Nxi2eta[i](xi, eta); }
    Type Nxieta2(const size_t i, const Type xi, const Type eta) const { return basic_Nxieta2[i](xi, eta); }
    Type Neta3  (const size_t i, const Type xi, const Type eta) const { return basic_Neta3  [i](xi, eta); }
*/

/*
        basic_Nxi2 = { [](const Type xi, const Type eta) { return 9./16. * (1.-eta) * (1. - 3.*xi + (1.+2.*p)*(-1.-eta)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.-eta) * (1. - eta + 2.*p*(-1.-eta) - 18.*xi); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.-eta) * (1. - eta + 2.*p*(-1.-eta) + 18.*xi); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.-eta) * (1. + 3.*xi + (1.+2.*p)*(-1.-eta));; },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+eta*eta); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+eta*eta); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.+eta) * (1. + 3.*xi + (1.+2.*p)*(-1.+eta));; },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+eta) * (1. + eta + 2.*p*(-1.+eta) + 18.*xi); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+eta) * (1. + eta + 2.*p*(-1.+eta) - 18.*xi); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.+eta) * (1. - 3.*xi + (1.+2.*p)*(-1.+eta));; },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+eta*eta); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+eta*eta); } },

        basic_Nxieta={ [](const Type xi, const Type eta) { return 1./32. * (-18.*eta + 36.*(1.+2.*p)*xi*eta - (18.*xi - (-10.+27.*eta*eta) - 27*xi*xi)); },
                      [](const Type xi, const Type eta) { return 9./32. * (2.*xi*(1.-eta*(1.+2.*p)) + 3. - 9.*xi*xi); },
                      [](const Type xi, const Type eta) { return 9./32. * (2.*xi*(1.-eta*(1.+2.*p)) - 3. + 9.*xi*xi); },
                      [](const Type xi, const Type eta) { return 1./32. * ( 18.*eta + 36.*(1.+2.*p)*xi*eta - (18.*xi + (-10.+27.*eta*eta) + 27*xi*xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * ( 3. - 9.*eta*eta + 2.*eta*(1. + xi*(1+2.*p))); },
                      [](const Type xi, const Type eta) { return -9./32. * (-3. + 9.*eta*eta + 2.*eta*(1. + xi*(1+2.*p))); },
                      [](const Type xi, const Type eta) { return 1./32. * ( 18.*eta + 36.*(1.+2.*p)*xi*eta + (18.*xi + (-10.+27.*eta*eta) + 27*xi*xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (2.*xi*(1.+eta*(1.+2.*p)) - 3. + 9.*xi*xi); },
                      [](const Type xi, const Type eta) { return -9./32. * (2.*xi*(1.+eta*(1.+2.*p)) + 3. - 9.*xi*xi); },
                      [](const Type xi, const Type eta) { return 1./32. * (-18.*eta + 36.*(1.+2.*p)*xi*eta + (18.*xi - (-10.+27.*eta*eta) - 27*xi*xi)); },
                      [](const Type xi, const Type eta) { return  9./32. * (-3. + 9.*eta*eta + 2.*eta*(1. - xi*(1+2.*p))); },
                      [](const Type xi, const Type eta) { return  9./32. * ( 3. - 9.*eta*eta + 2.*eta*(1. - xi*(1+2.*p))); } },

        basic_Neta2= { [](const Type xi, const Type eta) { return 9./16. * (1.-xi) * (1. - 3.*eta + (1.+2.*p) * (-1.-xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+xi*xi); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+xi*xi); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.+xi) * (1. - 3.*eta + (1.+2.*p) * (-1.+xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+xi) * (1. - 18.*eta + xi + 2.*p*(-1.+xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+xi) * (1. + 18.*eta + xi + 2.*p*(-1.+xi)); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.+xi) * (1. + 3.*eta + (1.+2.*p) * (-1.+xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+xi*xi); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.+2.*p) * (-1.+xi*xi); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.-xi) * (1. + 3.*eta + (1.+2.*p) * (-1.-xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.-xi) * (1. + 18.*eta - xi + 2.*p*(-1.-xi)); },
                      [](const Type xi, const Type eta) { return -9./32. * (1.-xi) * (1. - 18.*eta - xi + 2.*p*(-1.-xi)); } },

        basic_Nxi3 = { [](const Type xi, const Type eta) { return -27./16. * (1-eta); },
                      [](const Type xi, const Type eta) { return 81./16. * (1.-eta); },
                      [](const Type xi, const Type eta) { return -81./16. * (1.-eta); },
                      [](const Type xi, const Type eta) { return  27./16. * (1-eta); },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return  27./16. * (1+eta); },
                      [](const Type xi, const Type eta) { return -81./16. * (1.+eta); },
                      [](const Type xi, const Type eta) { return 81./16. * (1.+eta); },
                      [](const Type xi, const Type eta) { return -27./16. * (1+eta); },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return 0.; } },

     basic_Nxi2eta = { [](const Type xi, const Type eta) { return -9./16. * (1. - 3.*xi - 2.*(1.+2.*p)*eta); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.-eta*(1.+2.*p) - 9.*xi); },
                      [](const Type xi, const Type eta) { return 9./16. * (1.-eta*(1.+2.*p) + 9.*xi); },
                      [](const Type xi, const Type eta) { return -9./16. * (1. + 3.*xi - 2.*(1.+2.*p)*eta); },
                      [](const Type xi, const Type eta) { return -9./16. * eta * (1.+2.*p); },
                      [](const Type xi, const Type eta) { return -9./16. * eta * (1.+2.*p); },
                      [](const Type xi, const Type eta) { return  9./16. * (1. + 3.*xi + 2.*(1.+2.*p)*eta); },
                      [](const Type xi, const Type eta) { return -9./16. * (1.+eta*(1.+2.*p) + 9.*xi); },
                      [](const Type xi, const Type eta) { return -9./16. * (1.+eta*(1.+2.*p) - 9.*xi); },
                      [](const Type xi, const Type eta) { return  9./16. * (1. - 3.*xi + 2.*(1.+2.*p)*eta); },
                      [](const Type xi, const Type eta) { return -9./16. * eta * (1.+2.*p); },
                      [](const Type xi, const Type eta) { return -9./16. * eta * (1.+2.*p); } },

       basic_Nxieta2={ [](const Type xi, const Type eta) { return -9./16. * (1. - 3.*eta - 2.*(1.+2.*p)*xi); },
                      [](const Type xi, const Type eta) { return -9./16. * xi * (1. + 2.*p); },
                      [](const Type xi, const Type eta) { return -9./16. * xi * (1. + 2.*p); },
                      [](const Type xi, const Type eta) { return  9./16. * (1. - 3.*eta + 2.*(1.+2.*p)*xi); },
                      [](const Type xi, const Type eta) { return -9./16. * (1. - 9.*eta + xi*(1.+2.*p)); },
                      [](const Type xi, const Type eta) { return -9./16. * (1. + 9.*eta + xi*(1.+2.*p)); },
                      [](const Type xi, const Type eta) { return  9./16. * (1. + 3.*eta + 2.*(1.+2.*p)*xi); },
                      [](const Type xi, const Type eta) { return -9./16. * xi * (1. + 2.*p); },
                      [](const Type xi, const Type eta) { return -9./16. * xi * (1. + 2.*p); },
                      [](const Type xi, const Type eta) { return -9./16. * (1. + 3.*eta - 2.*(1.+2.*p)*xi); },
                      [](const Type xi, const Type eta) { return  9./16. * (1. + 9.*eta - xi*(1.+2.*p)); },
                      [](const Type xi, const Type eta) { return  9./16. * (1. - 9.*eta - xi*(1.+2.*p)); } },

        basic_Neta3= { [](const Type xi, const Type eta) { return -27./16. * (1.-xi); },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return -27./16. * (1.+xi); },
                      [](const Type xi, const Type eta) { return  81./16. * (1.+xi); },
                      [](const Type xi, const Type eta) { return -81./16. * (1.+xi); },
                      [](const Type xi, const Type eta) { return  27./16. * (1.+xi); },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return 0.; },
                      [](const Type xi, const Type eta) { return  27./16. * (1.-xi); },
                      [](const Type xi, const Type eta) { return -81./16. * (1.-xi); },
                      [](const Type xi, const Type eta) { return  81./16. * (1.-xi);; } };
*/

/*
// Частичная специализация класса описанного выше, на случай кубических серендиповых элементов. Необходима для экспериментальных целей.
template<class Type>
class element_2d_integrate<Type, qubic_serendip> : public element_2d_integrate_base<Type>,
                                                   public element_2d<Type, qubic_serendip> {
    std::vector<Type> weights, N_in_quad_nodes, Nxi_in_quad_nodes, Neta_in_quad_nodes,
                      Nxi2InQuad, NxietaInQuad, Neta2InQuad, Nxi3InQuad, Nxi2etaInQuad, Nxieta2InQuad, Neta3InQuad;

public:
    element_2d_integrate(const quadrature_base<Type> &quadrature_xi, const quadrature_base<Type> &quadrature_eta) { set(quadrature_xi, quadrature_eta); }

    // Узлы в которых происходит расчёт получаются путём декартова произведения квадратур с учётом того,
    // что геометрия описывается таким образом, что параметрическую зависимость имеет верхняя и нижняя границы,
    // а правая и левая заданы константами.
    void set(const quadrature_base<Type> &quadrature_xi, const quadrature_base<Type> &quadrature_eta) override
    {
        Type jacobian_xi = (this-> boundary(side_2d::RIGHT, 0) - this ->boundary(side_2d::LEFT, 0)) /
                          (quadrature_xi.boundary(side_1d::RIGHT)    - quadrature_xi.boundary(side_1d::LEFT));
        std::vector<Type> jacobian_eta(quadrature_xi. nodes_count()),
                          xi(quadrature_xi. nodes_count()),
                          eta(quadrature_eta.nodes_count());

        weights.resize(quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        for(size_t i = 0; i < quadrature_xi.nodes_count(); ++i) {
            xi[i] = this->boundary(side_2d::LEFT, 0) + (quadrature_xi.node(i) - quadrature_xi.boundary(side_1d::LEFT)) * jacobian_xi;
            jacobian_eta[i] = (this  ->boundary(side_2d::UP, xi[i]) - this  ->boundary(side_2d::DOWN, xi[i])) /
                             (quadrature_eta.boundary(side_1d::RIGHT)     - quadrature_eta.boundary(side_1d::LEFT));
            for(size_t j = 0; j < quadrature_eta.nodes_count(); ++j) {
                eta[j] = this->boundary(side_2d::DOWN, xi[i]) + (quadrature_eta.node(j)-quadrature_eta.boundary(side_1d::LEFT)) * jacobian_eta[i];
                weights[i*quadrature_eta.nodes_count() + j] = quadrature_xi.weight(i) * jacobian_xi * quadrature_eta.weight(j) * jacobian_eta[i];
            }
        }

        N_in_quad_nodes.resize(int_base::qnodes_count() * this->nodes_count());
        Nxi_in_quad_nodes.resize(int_base::qnodes_count() * this->nodes_count());
        Neta_in_quad_nodes.resize(int_base::qnodes_count() * this->nodes_count());

        Nxi2InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        NxietaInQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Neta2InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Nxi3InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Nxi2etaInQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Nxieta2InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Neta3InQuad.resize(int_base::qnodes_count() * this->nodes_count());

        for(size_t i = 0; i < this->nodes_count(); ++i)
            for(size_t j = 0; j < quadrature_xi.nodes_count(); ++j)
                for(size_t k = 0; k < quadrature_eta.nodes_count(); ++k) {
                    N_in_quad_nodes   [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->N   (i, xi[j], eta[k]);
                    Nxi_in_quad_nodes [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi (i, xi[j], eta[k]);
                    Neta_in_quad_nodes[i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Neta(i, xi[j], eta[k]);

                    Nxi2InQuad   [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi2   (i, xi[j], eta[k]);
                    NxietaInQuad [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxieta (i, xi[j], eta[k]);
                    Neta2InQuad  [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Neta2  (i, xi[j], eta[k]);
                    Nxi3InQuad   [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi3   (i, xi[j], eta[k]);
                    Nxi2etaInQuad[i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi2eta(i, xi[j], eta[k]);
                    Nxieta2InQuad[i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxieta2(i, xi[j], eta[k]);
                    Neta3InQuad  [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Neta3  (i, xi[j], eta[k]);
                }
    }

    size_t int_base::qnodes_count() const override { return weights.size(); }

    Type weight(const size_t q) const override { return weights[q]; }

    Type qN   (const size_t i, const size_t q) const override { return N_in_quad_nodes   [i*int_base::qnodes_count() + q]; }
    Type qNxi (const size_t i, const size_t q) const override { return Nxi_in_quad_nodes [i*int_base::qnodes_count() + q]; }
    Type qNeta(const size_t i, const size_t q) const override { return Neta_in_quad_nodes[i*int_base::qnodes_count() + q]; }

    Type qNxi2   (const size_t i, const size_t q) const { return Nxi2InQuad   [i*int_base::qnodes_count() + q]; }
    Type qNxieta (const size_t i, const size_t q) const { return NxietaInQuad [i*int_base::qnodes_count() + q]; }
    Type qNeta2  (const size_t i, const size_t q) const { return Neta2InQuad  [i*int_base::qnodes_count() + q]; }
    Type qNxi3   (const size_t i, const size_t q) const { return Nxi3InQuad   [i*int_base::qnodes_count() + q]; }
    Type qNxi2eta(const size_t i, const size_t q) const { return Nxi2etaInQuad[i*int_base::qnodes_count() + q]; }
    Type qNxieta2(const size_t i, const size_t q) const { return Nxieta2InQuad[i*int_base::qnodes_count() + q]; }
    Type qNeta3  (const size_t i, const size_t q) const { return Neta3InQuad  [i*int_base::qnodes_count() + q]; }

    virtual ~element_2d_integrate() = default;
};
*/