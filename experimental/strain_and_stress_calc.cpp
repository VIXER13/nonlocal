// Старый вариант вычисления деформаций и напряжений
// После дополнительных верификаций будет удалено

// Получение деформаций в узлах сетки, путём их переинтерполяции из квадратурных узлов
// Данный кусок взят из другой программы. Я не до конца понимаю как это работает, на мой взгляд это работать не должно.
// Пока что будем считать, что сетка однородная и состоит из билинейных элементов. В будущем, я надеюсь, это будет исправлено.
std::array<std::vector<double>, 3> strains_calc(const mesh_2d<double> &mesh, const Eigen::VectorXd &u)
{
    std::vector<double> eps11(mesh.nodes_count()),
                        eps22(mesh.nodes_count()),
                        eps12(mesh.nodes_count());

    std::vector<uint8_t> repeating(mesh.nodes_count(), 0);
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ++repeating[mesh.node_number(el, i)];
    }

    Eigen::MatrixXd NQP(e->nodes_count(), e->nodes_count());
    for(size_t q = 0; q < e->nodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
            NQP(q, i) = e->qN(i, q);
    Eigen::MatrixXd NQPI = NQP.inverse();

    Eigen::VectorXd EpsiXX_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                EpsiYY_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                EpsiXY_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                EpsiXX_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),
	                EpsiYY_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),	
	                EpsiXY_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),
	                Ue             = Eigen::VectorXd::Zero(2 * e->nodes_count()),
	                Epsi           = Eigen::VectorXd::Zero(3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * e->nodes_count()),
                    ElementNodesCoord(e->nodes_count(), 2),
                    Ndx(2, e->nodes_count());
    Eigen::Matrix2d Jmatr;
    Eigen::RowVectorXi ElementNodesNumbers(e->nodes_count());
    std::vector<Eigen::MatrixXd> NGradArr(e->nodes_count(), Eigen::MatrixXd::Zero(2, e->nodes_count()));

    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        B = Eigen::MatrixXd::Zero(3, 2 * e->nodes_count());
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ElementNodesNumbers[i] = mesh.node_number(el, i);

        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            ElementNodesCoord(i, 0) = mesh.coord(ElementNodesNumbers[i], 0);
            ElementNodesCoord(i, 1) = mesh.coord(ElementNodesNumbers[i], 1);
            Ue(i * 2)     = u(ElementNodesNumbers(i) * 2);
			Ue(i * 2 + 1) = u(ElementNodesNumbers(i) * 2 + 1);
        }

        for(size_t q = 0; q < e->nodes_count(); ++q)
        {
            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                NGradArr[q](0, j) = e->qNxi(j, q);
                NGradArr[q](1, j) = e->qNeta(j, q);
            }

            Jmatr = NGradArr[q] * ElementNodesCoord;
			Ndx = Jmatr.inverse() * NGradArr[q];
			for(size_t k = 0; k < e->nodes_count(); ++k)
			{
				B(0, k * 2)     = Ndx(0, k);
				B(1, k * 2 + 1) = Ndx(1, k);
				B(2, k * 2)     = Ndx(1, k);
				B(2, k * 2 + 1) = Ndx(0, k);
			}
			Epsi =  B * Ue;
			EpsiXX_QP(q) = Epsi(0);
			EpsiYY_QP(q) = Epsi(1);
			EpsiXY_QP(q) = Epsi(2);
        }

        EpsiXX_Element = NQPI * EpsiXX_QP;
		EpsiYY_Element = NQPI * EpsiYY_QP;
		EpsiXY_Element = 0.5 * NQPI * EpsiXY_QP;

        for(size_t i = 0; i < e->nodes_count(); ++i)
		{
			eps11[ElementNodesNumbers(i)] += EpsiXX_Element(i);
			eps22[ElementNodesNumbers(i)] += EpsiYY_Element(i);
			eps12[ElementNodesNumbers(i)] += EpsiXY_Element(i);
		}	
    }

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        eps11[i] /= repeating[i];
        eps22[i] /= repeating[i];
        eps12[i] /= repeating[i];
    }

    return {std::move(eps11), std::move(eps22), std::move(eps12)};
}

// И этой функции я тоже не доверяю, потому что взял её из другой программы
std::array<std::vector<double>, 3> stress_calc(const mesh_2d<double> &mesh, const Eigen::VectorXd &u, const parameters<double> &params)
{
    std::vector<double> sigma11(mesh.nodes_count()),
                        sigma22(mesh.nodes_count()),
                        sigma12(mesh.nodes_count());

    std::vector<uint8_t> repeating(mesh.nodes_count(), 0);
    const finite_element::element_2d_integrate_base<double> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ++repeating[mesh.node_number(el, i)];
    }

    Eigen::MatrixXd NQP(e->nodes_count(), e->nodes_count());
    for(size_t q = 0; q < e->nodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
            NQP(q, i) = e->qN(i, q);
    Eigen::MatrixXd NQPI = NQP.inverse();

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
    D(0, 0) = D(1, 1) = params.E / (1. - params.nu*params.nu);
    D(0, 1) = D(1, 0) = params.nu * params.E / (1. - params.nu*params.nu);
    D(2, 2) = 0.5 * params.E / (1. + params.nu);

    Eigen::VectorXd SigmaXX_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaYY_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaXY_Element = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaXX_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),
	                SigmaYY_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),	
	                SigmaXY_QP      = Eigen::VectorXd::Zero(    e->nodes_count()),
	                Ue              = Eigen::VectorXd::Zero(2 * e->nodes_count()),
	                Epsi            = Eigen::VectorXd::Zero(3),
                    Sigma           = Eigen::VectorXd::Zero(3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2 * e->nodes_count()),
                    ElementNodesCoord(e->nodes_count(), 2),
                    Ndx(2, e->nodes_count());
    Eigen::Matrix2d Jmatr;
    Eigen::RowVectorXi ElementNodesNumbers(e->nodes_count());
    std::vector<Eigen::MatrixXd> NGradArr(e->nodes_count(), Eigen::MatrixXd::Zero(2, e->nodes_count()));

    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        B = Eigen::MatrixXd::Zero(3, 2 * e->nodes_count());
        for(size_t i = 0; i < e->nodes_count(); ++i)
            ElementNodesNumbers[i] = mesh.node_number(el, i);

        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            ElementNodesCoord(i, 0) = mesh.coord(ElementNodesNumbers[i], 0);
            ElementNodesCoord(i, 1) = mesh.coord(ElementNodesNumbers[i], 1);
            Ue(i * 2)     = u(ElementNodesNumbers(i) * 2);
			Ue(i * 2 + 1) = u(ElementNodesNumbers(i) * 2 + 1);
        }

        for(size_t q = 0; q < e->nodes_count(); ++q)
        {
            for(size_t j = 0; j < e->nodes_count(); ++j)
            {
                NGradArr[q](0, j) = e->qNxi(j, q);
                NGradArr[q](1, j) = e->qNeta(j, q);
            }

            Jmatr = NGradArr[q] * ElementNodesCoord;
			Ndx = Jmatr.inverse() * NGradArr[q];
			for(size_t k = 0; k < e->nodes_count(); ++k)
			{
				B(0, k * 2)     = Ndx(0, k);
				B(1, k * 2 + 1) = Ndx(1, k);
				B(2, k * 2)     = Ndx(1, k);
                B(2, k * 2 + 1) = Ndx(0, k);
			}

            Epsi =  B * Ue;
			Epsi(2) *= 0.5;
			Sigma = D * Epsi;
			
			SigmaXX_QP(q) = Sigma(0);
			SigmaYY_QP(q) = Sigma(1);
			SigmaXY_QP(q) = Sigma(2);
        }

        SigmaXX_Element = NQPI * SigmaXX_QP;
		SigmaYY_Element = NQPI * SigmaYY_QP;
		SigmaXY_Element = NQPI * SigmaXY_QP;

        for(size_t i = 0; i < e->nodes_count(); ++i)
		{
			sigma11[ElementNodesNumbers(i)] += SigmaXX_Element(i);
			sigma22[ElementNodesNumbers(i)] += SigmaYY_Element(i);
			sigma12[ElementNodesNumbers(i)] += SigmaXY_Element(i);
		}
    }

    for(size_t i = 0; i < mesh.nodes_count(); ++i)
    {
        sigma11[i] /= repeating[i];
        sigma22[i] /= repeating[i];
        sigma12[i] /= repeating[i];
    }

    return {std::move(sigma11), std::move(sigma22), std::move(sigma12)};
}