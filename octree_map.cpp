#include "octree_map.h"

using namespace std;

namespace octree {


constexpr int kMinPointNumInOneCude = 2;

MapPoint::MapPoint(const float x, const float y, const float z) {
    p_ << x, y, z;
}

MapPoint::MapPoint(const Eigen::Vector3f &p) : p_(p){}


Cube::Cube(const float sideLen, const Eigen::Vector3f &center)
    : sideLen_(sideLen)
    , center_(center) {}

void Cube::CalculateChildNodeParameters(float &halfLen, float &quarterLen, float &x1, float &x2,
            float &y1, float &y2, float &z1, float &z2) {
    halfLen = sideLen_ * 0.5;
    quarterLen = sideLen_ * 0.25;
    x1 = center_.x() - quarterLen;
    x2 = center_.x() + quarterLen;
    y1 = center_.y() - quarterLen;
    y2 = center_.y() + quarterLen;
    z1 = center_.z() - quarterLen;
    z2 = center_.z() + quarterLen;        
}

void Cube::InitCubesFull() {
    cube.resize(8);
    float halfLen, quarterLen, x1, x2, y1, y2, z1, z2;
    CalculateChildNodeParameters(halfLen, quarterLen, x1, x2, y1, y2, z1, z2);
    /****** 使用如下八地图点存储顺序，以加快坐标所属node索引，单层最多比较3次即可 ******
    * 只有叶子结点的node才是最终分辨率，用于存储相应的map point
    * x1 x1 x1 x1 | x2 x2 x2 x2
    * y1 y1 y2 y2 | y1 y1 y2 y2
    * z1 z2 z1 z2 | z1 z2 z1 z2
    */
    cube[0] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y1, z1});
    cube[1] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y1, z2});
    cube[2] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y2, z1});
    cube[3] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y2, z2});

    cube[4] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y1, z1});
    cube[5] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y1, z2});
    cube[6] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y2, z1});
    cube[7] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y2, z2});
}

void Cube::CreateOneCube(const int id) {
    if(cube.empty()) {
        cube.resize(8);
    }
    if(!cube[id]) {
        float halfLen, quarterLen, x1, x2, y1, y2, z1, z2;
        CalculateChildNodeParameters(halfLen, quarterLen, x1, x2, y1, y2, z1, z2);

        switch (id) {
            case 0:
                cube[0] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y1, z1});
                return;
            case 1:
                cube[1] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y1, z2});
                return;
            case 2:
                cube[2] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y2, z1});
                return;
            case 3:
                cube[3] = make_shared<Cube>(halfLen, Eigen::Vector3f{x1, y2, z2});
                return;
            case 4:
                cube[4] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y1, z1});
                return;
            case 5:
                cube[5] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y1, z2});
                return;
            case 6:
                cube[6] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y2, z1});
                return;
            case 7:
                cube[7] = make_shared<Cube>(halfLen, Eigen::Vector3f{x2, y2, z2});
                return;
            default:
                cerr << id  << " not in[0, 7], Error!" << endl;
                exit(-1);
        }
    }
}



Octree::Octree(const float sideLen, const float resolution, const Eigen::Vector3f &center, const std::string &fileName)
    : sideLen_(sideLen)
    , resolution_(resolution)
    {
        debugFile_.open(fileName, ios::out);
        root_ = make_shared<Cube>(sideLen, center);
    }

void Octree::CreateOctrteeFull() {
    const auto t1 = chrono::high_resolution_clock::now();
        CreateOctrteeFull(root_);
        const auto t2 = chrono::high_resolution_clock::now();
        chrono::duration<double> spendTime = t2 - t1;
        cout << fixed;
        cout << "create sideLen= " << sideLen_ << "m & resolution= " 
                  << resolution_ << "m & leaf node num= " << leafNum_ 
                  << " spend " << spendTime.count() << " seconds." << endl;
        cout << "log2(leafNum_)= " << log2(leafNum_) << " & pow(2, x)= " << pow(2, log2(leafNum_)) << endl;
}

// 递归构建八叉树
void Octree::CreateOctrteeFull(shared_ptr<Cube> node) {
    // 递归终止条件
    // 最终八叉树的分辨率一般不会严格等于预设分辨率
    if(node->sideLen_ <= resolution_) {
        ++leafNum_;
        debugFile_ << node->center_.x() << " " << node->center_.y() << " " << node->center_.z() << endl;
        return;
    }

    // 当前一步操作：初始化结点
    node->InitCubesFull();
    for(int i = 0; i < 8; ++i) {
        CreateOctrteeFull(node->cube[i]);
    }
}

int Octree::FindCurrentPointBelong2WhichChildNode(const Eigen::Vector3f &p, std::shared_ptr<Cube> node){
    const Eigen::Vector3f &ct = node->center_;
    int id = INT_MAX;
    if(p.x() <= ct.x()) {
        // 答案在 cubes[0]~cudes[3]
        if(p.y() <= ct.y()) {
            // 答案在 cubes[0]~cubes[1]
            if(p.z() <= ct.z()) {
                // 答案为 cubes[0]
                id = 0;
            } else {
                // 答案为 cubes[1]
                id = 1;
            }
        } else {
            // 答案在 cubes[2]~cubes[3]
            if(p.z() <= ct.z()) {
                // 答案为 cubes[2]
                id = 2;
            } else {
                // 答案为 cubes[3]
                id = 3;
            }
        }
    } else {
        // 答案在 cubes[4]~cubes[7]
        if(p.y() <= ct.y()) {
            // 答案在 cubes[4]~cubes[5]
            if(p.z() <= ct.z()) {
                // 答案为 cubes[4]
                id = 4;
            } else {
                // 答案为 cubes[5]
                id = 5;
            }
        } else {
            // 答案在 cubes[6]~cubes[7]
            if(p.z() <= ct.z()) {
                // 答案为 cubes[6]
                id = 6;
            } else {
                // 答案为 cubes[7]
                id = 7;
            }
        }
    }
    return id;
}


// 实现动态扩增的八叉树，只在有地图点的位置生成leaf，而不需要构建满八叉树
bool Octree::AddOnePoint(const Eigen::Vector3f &p) {
    return AddOnePoint(p, root_);
}

bool Octree::AddOnePoint(const Eigen::Vector3f &p, shared_ptr<Cube> node) {
    // 递归终止条件
    if(node->sideLen_ <= resolution_) {
        node->points.push_back(make_shared<MapPoint>(p));
        return true;
    }

    // =========================================================
    // 当前步骤： 找到具体需要初始化哪个node
    const int id = FindCurrentPointBelong2WhichChildNode(p, node);
    // 当第id个结点是空时，我们需要创建它
    node->CreateOneCube(id);

    // ========================================================
    // 当已经创建好节点时，调用return并递归遍历
    return AddOnePoint(p, node->cube[id]);
}


shared_ptr<Cube> Octree::LeafPointBelong2(const MapPoint &p) {
    const Eigen::Vector3f &_p = p.p_;
    return LeafPointBelong2(_p, root_);
}

// 递归索引任意地图点所在的叶子结点
shared_ptr<Cube> Octree::LeafPointBelong2(const Eigen::Vector3f &p, shared_ptr<Cube> node){
    // 递归终止条件
    if(!node) {
        // 这是由于动态创建的octree可能会在相应索引位置不生成leaf node
        return node;
    }
    if(node->sideLen_ <= resolution_) {
        return node;
    }

    // ========================================================
    // 当前操作步骤
    const int id = FindCurrentPointBelong2WhichChildNode(p, node);

    const vector<shared_ptr<Cube>> &cubes = node->cube;
    if(cubes[id]) {
        cout << "see id: " << id << " | " << cubes[id]->center_.transpose() << " | sideLen: " 
             << cubes[id]->sideLen_ << endl;
    }

    //========================================================
    // 递归函数调用处return
    return LeafPointBelong2(p, cubes[id]);
}

void Octree::TraverseAllLeafNode(std::vector<std::shared_ptr<Cube> > &result, std::shared_ptr<Cube> node) {
    if(!node) {
        return;
    }
    if(node->sideLen_ <= resolution_) {
        result.push_back(node);
        return;
    }

    for(int i = 0; i < 8; ++i) {
        TraverseAllLeafNode(result, node->cube[i]);
    }
}

void Octree::TraverseAllLeafNode(std::vector<std::shared_ptr<Cube> > &result) {
    TraverseAllLeafNode(result, root_);
}


Octree::~Octree() {
    debugFile_.close();
}

}

using namespace octree;
int main(int argc, char** argv) {
    if (argc != 6) {
        cerr << endl 
             << "Usage: ./octree_map side_length resolution [p.x p.y p.z]2find!" << endl;
        return 1;
    }

    const float len = stof(argv[1]);
    const float res = stof(argv[2]);
    const Eigen::Vector3f ct{0, 0, 0};
    // const string debugPointCloudFileName = "point_cloud_full_octree.csv";
    const string debugPointCloudFileName = "point_cloud_dynamic_create_octree.csv";
    Octree octree(len, res, ct, debugPointCloudFileName);
    
    bool createFullOctree = false;
    if(createFullOctree) {
        // 创建满八叉树
        octree.CreateOctrteeFull();
    } else {
        // 读入点云数据文件，逐个创建八叉树结点
        int count = 0;
        ifstream fin;
        const string pointCloud2MakeOctreeFile("./point_cloud_2_make_octree.csv");
        fin.open(pointCloud2MakeOctreeFile.c_str());
            string input;
        while (getline(fin, input)) {
            if (input.empty()) {
                continue;
            }
            stringstream ss(input);
            double x, y, z;
            ss >> x >> y >> z;
            octree.AddOnePoint({x, y, z});
            ++count;
        }
        cout << "Add to octree map point num: " << count << endl;
        std::vector<std::shared_ptr<Cube> > result;
        octree.TraverseAllLeafNode(result);
        cout << "map point num after octree filtered: " << result.size()
             << " & culling point num: " << (count - result.size()) << endl;
        for(int i = 0; i < result.size(); ++i) {
            shared_ptr<Cube> leafNode = result[i];
            if(leafNode->points.size() >= kMinPointNumInOneCude) {
                octree.debugFile_ << leafNode->center_.x() << " " << leafNode->center_.y() << " " << leafNode->center_.z() << endl; 
            }
        }
    }

    // 寻找地图点所在位置
    const float x = stof(argv[3]);
    const float y = stof(argv[4]);
    const float z = stof(argv[5]);
    MapPoint p(x, y, z);
    shared_ptr<Cube> leaf = octree.LeafPointBelong2(p);
    if(leaf!=nullptr) {
        const Eigen::Vector3f lp = leaf->center_;
        cout << "point[" << x << " " << y << " " << z << "] is belong to leaf["
             << lp.x() << " " << lp.y() << " " << lp.z() << "] & their dist= " 
             << (p.p_ - lp).norm() << "< " << res << endl;
    } else {
        cout << "point[" << x << " " << y << " " << z << "] isn't belong to any existing leaf node" << endl;
    }

    return 0;
}
