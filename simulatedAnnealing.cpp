#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>

#include "Eigen/Dense"

class SimulatedAnnealing{
private:
    std::string atom_name;
    double eps;
    double sigma;
public:
    SimulatedAnnealing(std::string atom_name,double eps,double sigma){
        this->eps = eps;
        this->sigma = sigma;
    }

    Eigen::VectorXi selectPositions(int num, int N){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0,num);
        Eigen::VectorXi randomIndices(N);
        for (int i=0;i<N;i++){
            randomIndices(i) = dis(gen);
        }
        return randomIndices;
    }

    Eigen::MatrixXd displaceMolecules(Eigen::MatrixXd &pos,Eigen::MatrixXi indices, double d){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(-d,d);
        Eigen::MatrixXd newPos = pos;
        for (int j = 0;j<indices.size();j++){
            newPos(j,0) = pos(j,0) + dis(gen);
            newPos(j,1) = pos(j,1) + dis(gen);
            newPos(j,2) = pos(j,2) + dis(gen);
        }
        return newPos;
    }

    double distance(Eigen::MatrixXd &pos, int idx1, int idx2){
        double r = 0.0;
        Eigen::VectorXd pos1(3);
        Eigen::VectorXd pos2(3);
        for (int i=0;i<3;i++){
            pos1(i) = pos(idx1,i);
            pos2(i) = pos(idx2,i);
        }
        r = std::sqrt((pos1(0)-pos2(0))*(pos1(0)-pos2(0))+
                    (pos1(1)-pos2(1))*(pos1(1)-pos2(1))+
                    (pos1(2)-pos2(2))*(pos1(2)-pos2(2)));
        return r;
    }

    double TotalEnergy(int N, Eigen::MatrixXd &pos){
        double totalEnergy = 0.0;
        // for (int i=0;i<N;i++){
        //     std::cout << pos(i,0) << " , " << pos(i,1) << " , " << pos(i,2) << std::endl;
        // }
        for (int j = 0; j < N-1; j++){
            for (int k = j+1; k < N; k++){
                // std::cout << "(" << j <<","<< k << ")"<<std::endl;
                double r = distance(pos,j,k);
                if (r <= 12.0){
                    totalEnergy += 4*this->eps*(std::pow(this->sigma/r,12) - std::pow(this->sigma/r,6));
                }
            }
        }
        return totalEnergy;
    }

    // Generate a random number in a specified range
    double randomInRange(double min, double max) {
        return min + (max - min) * (static_cast<double>(rand()) / RAND_MAX);
    }

    int simulatedAnnealingLoop(double temperature,Eigen::MatrixXd &pos, double alpha){
        double totalEnergy = TotalEnergy(pos.rows(), pos);
        double newEnergy = totalEnergy;
        double minEnergy = totalEnergy;
        int N_iter = 1000;
        int accepted_count = 0;
        for (int iter=0;iter < N_iter; iter++){
            Eigen::VectorXi selectedPos = selectPositions(pos.rows(), 1);
            Eigen::MatrixXd newPos = displaceMolecules(pos,selectedPos,1.0);
            newEnergy = TotalEnergy(newPos.rows(), newPos);
            double deltaEnergy = newEnergy - totalEnergy;
            std::cout << "deltaEnergy : " << deltaEnergy <<std::endl;
            if (deltaEnergy < 0.0 || exp(-deltaEnergy/temperature) > randomInRange(0,1)){
                pos = newPos;
                totalEnergy = newEnergy;
                minEnergy = newEnergy;
                accepted_count++;
            }
            temperature *= alpha;
        }
        return accepted_count;
    }
};

Eigen::MatrixXd readCSVToEigenMatrix(const std::string &filename){
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening file");
    }

    std::vector<std::vector<double>> data;
    std::string line;

    while (std::getline(file,line)){
        std::stringstream ss(line);
        std::string item;
        std::vector<double> row;

        while(std::getline(ss,item,',')){
            row.push_back(std::stod(item));
        }
        data.push_back(row);
    }
    file.close();

    int rows = data.size();
    int cols = data[0].size();

    Eigen::MatrixXd mat(rows,cols);

    for (int i=0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            mat(i, j) = data[i][j];
        }
    }

    return mat;
}

int main(){
    SimulatedAnnealing simulatedAnnealing("CO2",130.0,3);
    Eigen::MatrixXd pos; //(2,3)
    // pos(0,0) = 2;
    // pos(0,1) = 0;
    // pos(0,2) = 3.4;
    // pos(1,0) = 2.3;
    // pos(1,1) = 4.42;
    // pos(1,2) = 5.52;
    pos = readCSVToEigenMatrix("pos.csv");
    std::cout << "The number of molecules is " << pos.rows() << std::endl;
    // std::cout << "Row " << 0 << " values: ";
    // for (int col = 0; col < pos.cols(); ++col) {
    //     std::cout << pos(0, col) << " ";
    // }
    double totalEnergy = simulatedAnnealing.TotalEnergy(pos.rows(),pos);
    // double distance = simulatedAnnealing.distance(pos,0,1);
    // std::cout << "The distance between points 0 and 1 is " << distance << std::endl; 
    std::cout << "The total energy before displacement is " << totalEnergy << " K" << std::endl; 

    // Eigen::MatrixXi selectedIndices = simulatedAnnealing.selectPositions(317,10);
    // std::cout << "Enter max displacement : ";
    // double maxDis = 0.0;
    // std::cin >>  maxDis;
    // simulatedAnnealing.displaceMolecules(pos,selectedIndices,maxDis);

    int accepted_count = simulatedAnnealing.simulatedAnnealingLoop(500,pos,0.9);
    std::cout << "Number of accepted moves : " << accepted_count << std::endl;
    double totalEnergyNew = simulatedAnnealing.TotalEnergy(pos.rows(),pos);
    std::cout << "The total energy after displacement is " << totalEnergyNew << " K" << std::endl; 


}