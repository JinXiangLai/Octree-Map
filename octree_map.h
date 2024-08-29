#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string.h>
#include <vector>

#include <Eigen/Dense>

namespace octree {

struct MapPoint {
    MapPoint(const float x, const float y, const float z);
    MapPoint(const Eigen::Vector3f &p);
    Eigen::Vector3f p_{0, 0, 0};
};

struct Cube {
    Cube(const float sideLen, const Eigen::Vector3f &center);
    void InitCubesFull();
    void CalculateChildNodeParameters(float &halfLen, float &quarterLen, float &x1, float &x2,
            float &y1, float &y2, float &z1, float &z2);
    void CreateOneCube(const int id);

    Eigen::Vector3f center_{0, 0, 0};
    float sideLen_ = 0.;
    std::vector<std::shared_ptr<Cube>> cube;

    // only used by the leaf nodes
    std::vector<std::shared_ptr<MapPoint>> points;
};

class Octree {
public:
    Octree(const float sideLen = 1000, const float resolution = 0.05, 
            const Eigen::Vector3f &center = Eigen::Vector3f::Zero(), 
            const std::string &fileName = "point_cloud.csv");
    ~Octree();
    void CreateOctrteeFull();
    std::shared_ptr<Cube> LeafPointBelong2(const MapPoint &p);
    bool AddOnePoint(const Eigen::Vector3f &p);
    void TraverseAllLeafNode(std::vector<std::shared_ptr<Cube> > &result);

    std::ofstream debugFile_;
private:
    void CreateOctrteeFull(std::shared_ptr<Cube> node);
    std::shared_ptr<Cube> LeafPointBelong2(const Eigen::Vector3f &p, std::shared_ptr<Cube> node);
    bool AddOnePoint(const Eigen::Vector3f &p, std::shared_ptr<Cube> node);
    int FindCurrentPointBelong2WhichChildNode(const Eigen::Vector3f &p, std::shared_ptr<Cube> node);
    void TraverseAllLeafNode(std::vector<std::shared_ptr<Cube> > &result, std::shared_ptr<Cube> node);

    std::shared_ptr<Cube> root_;
    float sideLen_ = 1000;
    float resolution_ = 0.05;
    int leafNum_ = 0;
};

}
