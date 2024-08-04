#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "Eigen/Dense"

using Eigen::MatrixXd;

Eigen::MatrixXi readCSVToEigenMatrix(const std::string &filename){
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening file");
    }

    std::vector<std::vector<int>> data;
    std::string line;

    while (std::getline(file,line)){
        std::stringstream ss(line);
        std::string item;
        std::vector<int> row;

        while(std::getline(ss,item,',')){
            row.push_back(std::stoi(item));
        }
        data.push_back(row);
    }
    file.close();

    int rows = data.size();
    int cols = data[0].size();

    Eigen::MatrixXi mat(rows,cols);

    for (int i=0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            mat(i, j) = data[i][j];
        }
    }

    return mat;
}

int main(){
    // Comma initialization
    // Eigen::Matrix3f nums;
    // nums << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    // std::cout << nums << std::endl;

    // Eigen::MatrixXd newArr;
    // newArr = MatrixXd::Constant(10,10,0);
    // std::cout << newArr.size() << std::endl;

    // Eigen::Matrix <float, 4, 4> newArr2;
    // newArr2.setZero();
    // std::cout << newArr2.size() << std::endl;

    // Eigen::MatrixXf results(4,4);
    // results << 8, 12, 17, 22,29, 101, 20, 91, 36, 42, 55,88, 10, 37 , 832, 812;
    // std::cout << "This is the original matrix results:\n" << results << std::endl;

    // results.transposeInPlace();
    // std::cout << "This is the new matrix results:\n" << results << std::endl;
     
    // sum(), prod(), mean(), minCoeff(), maxCoeff(),trace()
    try{
        Eigen::MatrixXi mat = readCSVToEigenMatrix("data.csv");
        // Eigen::Matrix2d mat;
        // mat.setOnes();
        std::cout << mat << std::endl;

        std::cout << "Result for sum() function: " << mat.sum() << std::endl;
        std::cout << "Result for prod() function: " << mat.prod() << std::endl;
        std::cout << "Result for minCoeff() function: " << mat.minCoeff() << std::endl;
        std::cout << "Result for maxCoeff() function: " << mat.maxCoeff() << std::endl;
        std::cout << "Result for trace() function: " << mat.trace() << std::endl;
        std::cout << "Result for mean() function: " << mat.mean() << std::endl;
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}