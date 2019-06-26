#include "FlexibleSolver.hpp"
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <fstream>
#include <iostream>
int
main()
{
    namespace pt = boost::property_tree;
    pt::ptree prm;

    {
        std::ifstream file("options.json");
        pt::read_json(file, prm);
        pt::write_json(std::cout, prm);
    }
    std::cout << "Hello, World!";
    constexpr int bz = 3;
    Dune::FlexibleSolver<bz> solver(prm);
    Dune::FlexibleSolver<2> solver2(prm);
    Dune::FlexibleSolver<1> solver1(prm);
    std::string matrixfile("matrix_istl.txt");
    std::string rhsfile("rhs_istl.txt");
    std::vector<double> res(9);
    // double* result;
    solver.solve(res.data(), matrixfile, rhsfile);
    std::cout << "**********Result*************" << std::endl;
    for (auto x : res) {
        std::cout << x << std::endl;
    }
    std::cout << "*********************************" << std::endl;


    int rows, cols, entries;
    std::vector<int> i;
    std::vector<int> j;
    std::vector<double> val;
    {
        std::ifstream file("matrix_istl.txt");
        if (!file) {
            throw std::runtime_error("Could not read");
        }
        std::string line;
        std::getline(file, line);
        std::cout << line << std::endl;
        std::getline(file, line);
        std::cout << line << std::endl;
        file >> rows;
        file >> cols;
        file >> entries;
        i.resize(entries);
        j.resize(entries);
        val.resize(entries);
        for (int kk = 0; kk < entries; ++kk) {
            file >> i[kk];
            --i[kk];
            file >> j[kk];
            --j[kk];
            file >> val[kk];
        }
    }
    std::vector<double> rhs;
    {
        std::ifstream file("rhs_istl.txt");
        if (!file) {
            throw std::runtime_error("Could not read");
        }
        std::string line;
        std::getline(file, line);
        std::cout << line << std::endl;
        std::getline(file, line);
        std::cout << line << std::endl;
        file >> rows;
        int tcols;
        file >> tcols;
        rhs.resize(rows);
        std::cout << "******** new rhs **********" << std::endl;
        for (int kk = 0; kk < rows; ++kk) {
            file >> rhs[kk];
            std::cout << rhs[kk] << std::endl;
        }
    }
    Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>> matrix;
    makeMatrixMarket(matrix, i, j, val, rows, cols, entries);

    {
        std::ofstream file("matrix_new.txt");
        writeMatrixMarket(matrix, file);
    }
    double tol = 1e-2;
    int maxiter = 200;
    if (rhs.size() % 3 == 0) {
        pt::ptree out;
        solver.solve(res.data(), i, j, val, rows, rhs, tol, maxiter, out);
        std::cout << "**********Result new*************" << std::endl;
        for (auto x : res) {
            std::cout << x << std::endl;
        }
    }
    if (rhs.size() % 2 == 0) {
        std::cout << "*********************************" << std::endl;
        pt::ptree out;
        solver2.solve(res.data(), i, j, val, rows, rhs, tol, maxiter, out);
        std::cout << "**********Result new*************" << std::endl;
        for (auto x : res) {
            std::cout << x << std::endl;
        }
    }
    {
        pt::ptree out;
        std::cout << "*********************************" << std::endl;
        solver1.solve(res.data(), i, j, val, rows, rhs, tol, maxiter, out);
        std::cout << "**********Result new*************" << std::endl;
    }
    for (auto x : res) {
        std::cout << x << std::endl;
    }
    std::cout << "*********************************" << std::endl;

    return 0;
}
